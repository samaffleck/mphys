#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "mphys/mesh.hpp"
#include "mphys/models/spme.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/transient_solver.hpp"

namespace {

using mphys::models::ElectrolyteMesh;
using mphys::models::MakeSpmeElectrolyteMesh;
using mphys::models::SpmeModel;
using mphys::models::SpmeParameters;
using mphys::models::SpmParameters;

// Volume-weighted average of a particle field over the spherical mesh.
double SphericalAverage(const mphys::Field& c, const mphys::Mesh1D& mesh) {
  double num = 0.0, den = 0.0;
  for (int i = 0; i < mesh.n_cells; ++i) {
    const double rl = mesh.face_positions[i];
    const double rr = mesh.face_positions[i + 1];
    const double vol = (rr * rr * rr - rl * rl * rl) / 3.0;  // 4π/3 cancels
    num += vol * c[i];
    den += vol;
  }
  return num / den;
}

// Porosity-weighted electrolyte mass ∫ ε c_e dx (Cartesian volume = dx).
double ElectrolyteMass(const mphys::Field& ce, const ElectrolyteMesh& em,
                       const SpmeParameters& p) {
  double m = 0.0;
  for (int i = 0; i < em.mesh.n_cells; ++i) {
    const double eps = em.region[i] == 0 ? p.eps_e_n
                     : em.region[i] == 1 ? p.eps_e_s
                                         : p.eps_e_p;
    m += eps * ce[i] * em.mesh.dx[i];
  }
  return m;
}

// Drive an SPMe to t_end and return the final field set.
std::vector<mphys::Field> RunToEnd(SpmeModel& model, double t_end,
                                   double max_dt = 20.0) {
  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-3;
  opts.maximum_time_step = max_dt;
  opts.tolerance.relative = 1e-9;
  opts.tolerance.absolute = 1e-9;

  mphys::SunContext sunctx;
  mphys::TransientSolver solver(model, opts, sunctx);

  std::vector<mphys::Field> ff;
  solver.Solve(0.0, t_end,
               [&](double, const std::vector<mphys::Field>& f,
                   const std::vector<double>&) { ff = f; });
  return ff;
}

}  // namespace

// The particle equations are unchanged from the SPM, so lithium couliometry in
// the negative particle must still hold exactly: d⟨c_n⟩/dt = -I/(A F ε_n L_n).
TEST(Spme, ParticleCouliometry) {
  SpmeParameters p;          // LG M50 defaults
  p.core.I = 2.0;
  const double t_end = 1500.0;

  auto em = MakeSpmeElectrolyteMesh(p, 20, 12, 20);
  auto particle_mesh = mphys::MakeUniformMesh1D(
      0.0, 1.0, em.mesh.n_cells, mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(em.mesh.n_cells);
  SpmeModel model(particle_mesh, em, sv, p);

  const auto ff = RunToEnd(model, t_end);

  const double rate = -p.core.I /
                      (p.core.A * SpmParameters::F * p.core.eps_n * p.core.L_n);
  const double expected = p.core.cn0 + rate * t_end;
  const double avg = SphericalAverage(ff[model.c_n_index()], particle_mesh);
  EXPECT_NEAR(avg, expected, 1e-2);
}

// The electrolyte source integrates to zero and the boundaries are no-flux, so
// the porosity-weighted electrolyte mass ⟨ε c_e⟩ is conserved for all time.
TEST(Spme, ElectrolyteMassConservation) {
  SpmeParameters p;
  p.core.I = 3.0;
  const double t_end = 1000.0;

  auto em = MakeSpmeElectrolyteMesh(p, 25, 15, 25);
  auto particle_mesh = mphys::MakeUniformMesh1D(
      0.0, 1.0, em.mesh.n_cells, mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(em.mesh.n_cells);
  SpmeModel model(particle_mesh, em, sv, p);

  mphys::Field ce0("c_e", em.mesh.n_cells, p.ce0);
  const double mass0 = ElectrolyteMass(ce0, em, p);

  const auto ff = RunToEnd(model, t_end);
  const double mass1 = ElectrolyteMass(ff[model.c_e_index()], em, p);

  EXPECT_NEAR(mass1, mass0, 1e-6 * mass0);
}

// Under constant current the electrolyte relaxes to a genuine steady profile
// satisfying  d/dx(D_e^eff dc_e/dx) = -S(x)  with no-flux ends.  Integrating the
// flux N(x) = ∫_0^x S and then c_e' = -N/D_e^eff gives a piecewise quadratic /
// linear / quadratic profile, fixed by conservation of ⟨ε c_e⟩.  The numerical
// solution must converge to this analytic profile.
TEST(Spme, ElectrolyteSteadyProfile) {
  SpmeParameters p;
  p.core.I = 0.5;            // C/10-ish: keeps particles valid over a long run
  const double t_end = 10000.0;

  const int nn = 40, ns = 30, np = 40;
  auto em = MakeSpmeElectrolyteMesh(p, nn, ns, np);
  auto particle_mesh = mphys::MakeUniformMesh1D(
      0.0, 1.0, em.mesh.n_cells, mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(em.mesh.n_cells);
  SpmeModel model(particle_mesh, em, sv, p);

  const auto ff = RunToEnd(model, t_end, 50.0);
  const auto& ce = ff[model.c_e_index()];

  // Analytic steady profile (base solution with base(0)=0).
  const double Ln = p.core.L_n, Ls = p.L_s, Lp = p.core.L_p;
  const double i_app = p.core.i_app();
  const double s = (1.0 - p.t_plus) * i_app / SpmParameters::F;  // plateau flux N0
  const double De_n = p.De_n(), De_s = p.De_s(), De_p = p.De_p();

  const double bn = -s * Ln / (2.0 * De_n);          // base at x = Ln
  const double bs = bn - s * Ls / De_s;              // base at x = Ln+Ls
  auto base = [&](double x) {
    if (x <= Ln) {
      return -(s / Ln) * x * x / (2.0 * De_n);
    } else if (x <= Ln + Ls) {
      return bn - s * (x - Ln) / De_s;
    } else {
      const double xi = x - (Ln + Ls);
      const double Sp = -s / Lp;
      return bs - (s * xi + Sp * xi * xi / 2.0) / De_p;
    }
  };

  // Fix the additive constant via ⟨ε·c_e⟩ conservation (initial value p.ce0).
  double num = 0.0, den = 0.0;
  for (int i = 0; i < em.mesh.n_cells; ++i) {
    const double eps = em.region[i] == 0 ? p.eps_e_n
                     : em.region[i] == 1 ? p.eps_e_s
                                         : p.eps_e_p;
    num += eps * base(em.mesh.cell_centres[i]) * em.mesh.dx[i];
    den += eps * em.mesh.dx[i];
  }
  const double C = p.ce0 - num / den;

  double max_err = 0.0;
  for (int i = 0; i < em.mesh.n_cells; ++i) {
    const double analytic = base(em.mesh.cell_centres[i]) + C;
    max_err = std::max(max_err, std::abs(ce[i] - analytic));
  }
  EXPECT_LT(max_err, 1.0);   // mol/m³, vs a ~70 mol/m³ profile variation
}

// The reported voltage must equal the closed-form SPMe expression at the
// solved state, sit in a physical Li-ion range, and the electrolyte corrections
// must actually shift it relative to the bare SPM open-circuit/reaction voltage.
TEST(Spme, VoltageConsistencyAndCorrections) {
  SpmeParameters p;
  p.core.I = 5.0;
  const double t_end = 200.0;

  auto em = MakeSpmeElectrolyteMesh(p, 20, 12, 20);
  auto particle_mesh = mphys::MakeUniformMesh1D(
      0.0, 1.0, em.mesh.n_cells, mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(em.mesh.n_cells);
  SpmeModel model(particle_mesh, em, sv, p);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-3;
  opts.maximum_time_step = 5.0;
  opts.tolerance.relative = 1e-8;
  opts.tolerance.absolute = 1e-8;

  mphys::SunContext sunctx;
  mphys::TransientSolver solver(model, opts, sunctx);

  std::vector<mphys::Field> ff;
  std::vector<double> fa;
  solver.Solve(0.0, t_end,
               [&](double, const std::vector<mphys::Field>& f,
                   const std::vector<double>& a) { ff = f; fa = a; });

  const double gn = -p.core.R_n * p.core.R_n * p.core.i_app() /
                    (3.0 * p.core.D_n * SpmParameters::F * p.core.eps_n * p.core.L_n);
  const double gp = +p.core.R_p * p.core.R_p * p.core.i_app() /
                    (3.0 * p.core.D_p * SpmParameters::F * p.core.eps_p * p.core.L_p);
  const double cn_s = model.SurfaceValue(ff[model.c_n_index()], gn);
  const double cp_s = model.SurfaceValue(ff[model.c_p_index()], gp);
  const double ce_n = model.RegionAverage(ff[model.c_e_index()], 0);
  const double ce_p = model.RegionAverage(ff[model.c_e_index()], 2);

  const double v_expected =
      mphys::models::SpmeVoltage(p, cn_s, cp_s, ce_n, ce_p);
  const double v_solver = fa[model.voltage_index()];

  EXPECT_NEAR(v_solver, v_expected, 1e-6);
  EXPECT_GT(v_solver, 2.0);
  EXPECT_LT(v_solver, 4.5);

  // Removing all electrolyte corrections (η_c, ΔΦ_e, ΔΦ_s) yields a different
  // voltage, confirming the corrections are active.
  SpmeParameters p_no_elec = p;
  p_no_elec.t_plus = 1.0;          // kills η_c
  p_no_elec.kappa_e = 1e30;        // kills ΔΦ_e
  p_no_elec.sigma_n = 1e30;        // kills ΔΦ_s
  p_no_elec.sigma_p = 1e30;
  const double v_bare =
      mphys::models::SpmeVoltage(p_no_elec, cn_s, cp_s, ce_n, ce_p);
  EXPECT_GT(std::abs(v_solver - v_bare), 1e-3);
}
