#include <memory>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/models/darcy_packed_bed.hpp"
#include "mphys/models/model_registry.hpp"
#include "mphys/models/transport.hpp"

namespace mphys {
namespace {

// Standard Dirichlet/Neumann boundary: option 0 = Dirichlet, 1 = Neumann.
BoundaryCondition StdBc(const BcChoice& c) {
  return c.option == 0 ? DirichletBc(c.value) : NeumannBc(c.value);
}

BcSlot DirichletNeumannSlot(const char* key, const char* label,
                            int default_option, double default_value) {
  return BcSlot{key, label,
                {BcOption{"Dirichlet", ""}, BcOption{"Neumann", ""}},
                default_option, default_value};
}

// ── Transport :: Convection-Diffusion-Reaction ──────────────────────────────
ModelInfo ConvDiffReactionInfo() {
  ModelInfo m;
  m.id = "transport.conv_diff_reaction";
  m.package = "Transport";
  m.name = "Convection-Diffusion-Reaction";
  m.description = "Transient transport of a scalar with diffusion, "
                  "first-order upwind convection and linear reaction.";
  m.solver = SolverKind::kTransient;
  m.schema.params = {
      {"D", "Diffusivity", "m\xc2\xb2/s", 1.0, ParamScope::kPerDomain,
       FaceInterp::kHarmonic},
      {"u", "Velocity", "m/s", 0.0, ParamScope::kPerDomain,
       FaceInterp::kArithmetic},
      {"k", "Reaction rate", "1/s", 0.0, ParamScope::kPerDomain,
       FaceInterp::kNone},
  };
  m.schema.boundaries = {
      DirichletNeumannSlot("left", "Boundary 1", 0, 1.0),
      DirichletNeumannSlot("right", "Boundary 2", 1, 0.0),
  };
  m.factory = [](ModelBuildContext& ctx) -> std::unique_ptr<Model> {
    auto D_face = ctx.FaceParam("D", FaceInterp::kHarmonic);
    auto u_face = ctx.FaceParam("u", FaceInterp::kArithmetic);
    auto k_cell = ctx.CellParam("k");
    return std::make_unique<models::ConvDiffReactionModel>(
        ctx.mesh(), ctx.sv(), std::move(D_face), std::move(u_face),
        std::move(k_cell), StdBc(ctx.Bc("left")), StdBc(ctx.Bc("right")));
  };
  return m;
}

// ── Transport :: Steady Diffusion ───────────────────────────────────────────
ModelInfo SteadyDiffusionInfo() {
  ModelInfo m;
  m.id = "transport.steady_diffusion";
  m.package = "Transport";
  m.name = "Steady Diffusion";
  m.description = "Steady-state diffusion with spatially-varying diffusivity.";
  m.solver = SolverKind::kSteady;
  m.schema.params = {
      {"D", "Diffusivity", "m\xc2\xb2/s", 1.0, ParamScope::kPerDomain,
       FaceInterp::kHarmonic},
  };
  m.schema.boundaries = {
      DirichletNeumannSlot("left", "Boundary 1", 0, 1.0),
      DirichletNeumannSlot("right", "Boundary 2", 1, 0.0),
  };
  m.factory = [](ModelBuildContext& ctx) -> std::unique_ptr<Model> {
    auto D_face = ctx.FaceParam("D", FaceInterp::kHarmonic);
    return std::make_unique<models::SteadyDiffusionModel>(
        ctx.mesh(), ctx.sv(), std::move(D_face), StdBc(ctx.Bc("left")),
        StdBc(ctx.Bc("right")));
  };
  return m;
}

// ── Fluid Flow :: Darcy Packed Bed ──────────────────────────────────────────
ModelInfo DarcyPackedBedInfo() {
  ModelInfo m;
  m.id = "fluid_flow.darcy_packed_bed";
  m.package = "Fluid Flow";
  m.name = "Darcy Packed Bed";
  m.description = "Darcy's law coupled to an ideal-gas mass balance through a "
                  "porous packed bed.";
  m.solver = SolverKind::kTransient;
  m.schema.params = {
      {"kappa", "Permeability", "m\xc2\xb2", 1e-10, ParamScope::kPerDomain,
       FaceInterp::kNone},
      {"mu", "Dynamic viscosity", "Pa\xc2\xb7s", 1.8e-5, ParamScope::kPerDomain,
       FaceInterp::kNone},
  };
  // Each boundary specifies either a Pressure or a Velocity.
  BcSlot left{"left", "Boundary 1",
              {BcOption{"Pressure", "Pa"}, BcOption{"Velocity", "m/s"}}, 0, 2.0};
  BcSlot right{"right", "Boundary 2",
               {BcOption{"Pressure", "Pa"}, BcOption{"Velocity", "m/s"}}, 0, 1.0};
  m.schema.boundaries = {left, right};
  m.factory = [](ModelBuildContext& ctx) -> std::unique_ptr<Model> {
    const auto& mesh = ctx.mesh();
    auto kappa_cell = ctx.CellParam("kappa");
    auto mu_cell = ctx.CellParam("mu");

    int nd = ctx.n_domains();
    double kL = ctx.DomainParam(0, "kappa"), muL = ctx.DomainParam(0, "mu");
    double kR = ctx.DomainParam(nd - 1, "kappa");
    double muR = ctx.DomainParam(nd - 1, "mu");

    BcChoice left_bc = ctx.Bc("left");
    BcChoice right_bc = ctx.Bc("right");

    // Pressure BC (option 0) → Dirichlet on P; Velocity BC (option 1) → Neumann
    // on P (from Darcy's law).
    auto make_P_bc = [](const BcChoice& bc, double kp, double mp) {
      return bc.option == 0 ? DirichletBc(bc.value)
                            : NeumannBc(-mp * bc.value / kp);
    };
    // Velocity BC → Dirichlet on u; Pressure BC → Neumann (free) on u.
    auto make_u_bc = [](const BcChoice& bc) {
      return bc.option == 1 ? DirichletBc(bc.value) : NeumannBc(0.0);
    };

    auto P_lbc = make_P_bc(left_bc, kL, muL);
    auto P_rbc = make_P_bc(right_bc, kR, muR);
    auto u_lbc = make_u_bc(left_bc);
    auto u_rbc = make_u_bc(right_bc);

    // Initialise P with a linear profile so IDACalcIC has a consistent start.
    double P_scalar = (left_bc.option == 0)    ? left_bc.value
                      : (right_bc.option == 0) ? right_bc.value
                                               : 101325.0;
    double P_L_val = (left_bc.option == 0) ? left_bc.value : P_scalar;
    double P_R_val = (right_bc.option == 0) ? right_bc.value : P_scalar;
    std::vector<double> P_init(mesh.n_cells);
    double x_L = mesh.cell_centres.front(), x_R = mesh.cell_centres.back();
    double span = (x_R > x_L) ? (x_R - x_L) : 1.0;
    for (int i = 0; i < mesh.n_cells; ++i)
      P_init[i] = P_L_val + (P_R_val - P_L_val) *
                                (mesh.cell_centres[i] - x_L) / span;

    return std::make_unique<models::DarcyPackedBedModel>(
        mesh, ctx.sv(), std::move(kappa_cell), std::move(mu_cell), P_init,
        P_lbc, P_rbc, u_lbc, u_rbc);
  };
  return m;
}

// ── Lithium-Ion Battery :: Single Particle Model ────────────────────────────
// SPM builds its own spherical mesh, exposes many global parameters, and has a
// bespoke voltage-vs-time results view, so it opts into custom GUI rather than
// the generic schema-driven path. It is catalogued here only so it appears in
// the package tree; the GUI owns its construction and rendering.
ModelInfo SingleParticleInfo() {
  ModelInfo m;
  m.id = "liion.spm";
  m.package = "Lithium-Ion Battery";
  m.name = "Single Particle Model";
  m.description = "Two spherical electrode particles with Fick diffusion and "
                  "Butler-Volmer kinetics; terminal voltage output (PyBaMM SPM).";
  m.solver = SolverKind::kTransient;
  m.custom_gui = true;
  return m;
}

}  // namespace

void RegisterBuiltinModels(ModelRegistry& registry) {
  registry.Register(ConvDiffReactionInfo());
  registry.Register(SteadyDiffusionInfo());
  registry.Register(DarcyPackedBedInfo());
  registry.Register(SingleParticleInfo());
}

}  // namespace mphys
