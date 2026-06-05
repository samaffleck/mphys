#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/fvm_mesh.hpp"
#include "mphys/mesh_model.hpp"
#include "mphys/mesh_steady_solver.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/topology.hpp"

namespace {

constexpr double pi = std::numbers::pi;

// 2D Poisson: laplacian(u) = f on [0,1]^2, u = 0 on the boundary.
// Manufactured solution u_exact = sin(pi x) sin(pi y), so
//   f = laplacian(u_exact) = -2 pi^2 sin(pi x) sin(pi y).
class PoissonModel : public mphys::MeshModel {
 public:
  explicit PoissonModel(const mphys::Mesh& mesh) : MeshModel(mesh) {
    u_ = AddField("u", 0.0);
    SetBcs(u_, {mphys::DirichletBc(0.0), mphys::DirichletBc(0.0),
                mphys::DirichletBc(0.0), mphys::DirichletBc(0.0)});
    // Precompute the source at cell centres.
    f_.resize(mesh.NCells());
    for (int c = 0; c < mesh.NCells(); ++c) {
      const auto& x = mesh.cells[c].centroid;
      f_[c] = -2.0 * pi * pi * std::sin(pi * x[0]) * std::sin(pi * x[1]);
    }
  }

  void Residual(double /*t*/, const std::vector<std::vector<double>>& y,
                const std::vector<std::vector<double>>& /*ydot*/,
                std::vector<std::vector<double>>& rr) override {
    const auto lap = mphys::fvm::Laplacian(y[u_], 1.0, mesh_, bcs(u_));
    for (int c = 0; c < mesh_.NCells(); ++c) rr[u_][c] = lap[c] - f_[c];
  }

 private:
  int u_ = 0;
  std::vector<double> f_;
};

double SolveAndMeasureError(int n) {
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, n, 0.0, 1.0, n);
  PoissonModel model(mesh);

  mphys::SolverOptions opts;
  opts.tolerance.absolute = 1e-11;
  opts.tolerance.relative = 1e-12;

  mphys::SunContext sunctx;
  mphys::MeshSteadySolver solver(model, opts, sunctx);
  solver.Solve();

  const auto& u = model.fields()[0];
  double err = 0.0;
  for (int c = 0; c < mesh.NCells(); ++c) {
    const auto& x = mesh.cells[c].centroid;
    const double exact = std::sin(pi * x[0]) * std::sin(pi * x[1]);
    err = std::max(err, std::abs(u[c] - exact));
  }
  return err;
}

}  // namespace

// End-to-end: build a 2D Poisson model on the face-based mesh and solve it with
// the matrix-free Newton-Krylov solver. The recovered solution must match the
// manufactured solution to within the expected discretisation error.
TEST(MeshSolve, Poisson2DAccuracy) {
  const double err = SolveAndMeasureError(32);
  // 32x32 mesh: O(h^2) ~ (1/32)^2 ~ 1e-3; allow comfortable slack.
  EXPECT_LT(err, 5e-3) << "L_inf error = " << err;
}

// The end-to-end solve must converge at ~2nd order under mesh refinement.
TEST(MeshSolve, Poisson2DSecondOrder) {
  const double e1 = SolveAndMeasureError(16);
  const double e2 = SolveAndMeasureError(32);
  const double rate = std::log(e1 / e2) / std::log(2.0);
  EXPECT_GT(rate, 1.8) << "observed convergence rate " << rate
                       << " (e16=" << e1 << ", e32=" << e2 << ")";
}
