#include <print>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/fvm_operators.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/sim_result.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/transient_solver.hpp"

// 1D transient convection-diffusion-reaction in a packed bed reactor.
//
// PDE: dc/dt = -u * dc/dx + D * d²c/dx² - k * c
//
// Domain: x in [0, 1] m
// BCs:    c(0, t) = 1.0   (Dirichlet inlet)
//         dc/dx(1, t) = 0  (Neumann outlet)
// IC:     c(x, 0) = 0
class ReactorModel : public mphys::Model {
 public:
  ReactorModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv,
               double D, double u, double k)
      : Model(mesh, sv), D_(D), u_(u), k_(k) {
    c_ = AddField("c", 0.0);
    SetBcs(c_, {mphys::DirichletBc(1.0), mphys::NeumannBc(0.0)});
  }

  void Residual(double /*t*/,
                const std::vector<mphys::Field>& y,
                const std::vector<mphys::Field>& ydot,
                const std::vector<double>& /*y_alg*/,
                const std::vector<double>& /*ydot_alg*/,
                std::vector<mphys::Field>& rr,
                std::vector<double>& /*rr_alg*/) override {
    rr[c_] = mphys::fvm::Ddt(ydot[c_])
           + mphys::fvm::Convection(y[c_], u_, mesh_, bcs_[c_])
           - mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_])
           + y[c_] * k_;
  }

 private:
  int c_ = 0;  // concentration field index
  double D_;   // diffusivity [m²/s]
  double u_;   // velocity [m/s]
  double k_;   // reaction rate [1/s]
};

int main() {
  constexpr int    kNCells = 100;
  constexpr double kD = 1e-4;
  constexpr double kU = 0.01;
  constexpr double kK = 0.1;
  constexpr double kTEnd = 50.0;

  auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, kNCells);

  mphys::StateVector sv(mesh.n_cells);
  ReactorModel model(mesh, sv, kD, kU, kK);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-3;
  opts.maximum_time_step = 1.0;
  opts.tolerance.relative = 1e-6;
  opts.tolerance.absolute = 1e-8;

  mphys::SunContext sunctx;
  mphys::SimResult result;

  mphys::TransientSolver solver(model, opts, sunctx);
  solver.Solve(0.0, kTEnd,
               [&](double t, const std::vector<mphys::Field>& fields,
                   const std::vector<double>& alg) {
                 result.Record(t, fields, alg);
               });

  std::println("Simulation complete. {} snapshots recorded.", result.snapshots.size());

  // Print the final concentration profile at a few points.
  const auto& final_c = result.snapshots.back().fields[0];
  std::println("Final concentration profile (every 10th cell):");
  for (int i = 0; i < mesh.n_cells; i += 10) {
    std::println("  x = {:.3f}  c = {:.6f}", mesh.cell_centres[i], final_c[i]);
  }

  return 0;
}
