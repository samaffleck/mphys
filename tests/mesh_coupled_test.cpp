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
constexpr double kAlpha = 5.0;  // cross-coupling strength

// Two linearly-coupled reaction-diffusion fields solved together in a single
// MeshModel:
//   laplacian(u) = alpha * v + g_u
//   laplacian(v) = alpha * u + g_v
// Manufactured solution (both zero on the boundary):
//   u_exact = sin(pi x) sin(pi y)
//   v_exact = x(1-x) y(1-y)
// The sources g_u, g_v are chosen so the exact pair satisfies the system, which
// lets us assert the coupled solve converges to it at 2nd order.
double UExact(double x, double y) { return std::sin(pi * x) * std::sin(pi * y); }
double VExact(double x, double y) { return x * (1.0 - x) * y * (1.0 - y); }
double LapUExact(double x, double y) { return -2.0 * pi * pi * UExact(x, y); }
double LapVExact(double x, double y) {
  return -2.0 * (y * (1.0 - y) + x * (1.0 - x));
}

class CoupledModel : public mphys::MeshModel {
 public:
  explicit CoupledModel(const mphys::Mesh& mesh) : MeshModel(mesh) {
    u_ = AddField("u", 0.0);
    v_ = AddField("v", 0.0);
    const std::vector<mphys::PatchBc> zero_dirichlet(
        mesh.patches.size(), mphys::DirichletBc(0.0));
    SetBcs(u_, zero_dirichlet);
    SetBcs(v_, zero_dirichlet);

    g_u_.resize(mesh.NCells());
    g_v_.resize(mesh.NCells());
    for (int c = 0; c < mesh.NCells(); ++c) {
      const auto& p = mesh.cells[c].centroid;
      g_u_[c] = LapUExact(p[0], p[1]) - kAlpha * VExact(p[0], p[1]);
      g_v_[c] = LapVExact(p[0], p[1]) - kAlpha * UExact(p[0], p[1]);
    }
  }

  void Residual(const std::vector<std::vector<double>>& y,
                std::vector<std::vector<double>>& rr) override {
    const auto lap_u = mphys::fvm::Laplacian(y[u_], 1.0, mesh_, bcs(u_));
    const auto lap_v = mphys::fvm::Laplacian(y[v_], 1.0, mesh_, bcs(v_));
    for (int c = 0; c < mesh_.NCells(); ++c) {
      rr[u_][c] = lap_u[c] - kAlpha * y[v_][c] - g_u_[c];
      rr[v_][c] = lap_v[c] - kAlpha * y[u_][c] - g_v_[c];
    }
  }

  int u() const { return u_; }
  int v() const { return v_; }

 private:
  int u_ = 0, v_ = 1;
  std::vector<double> g_u_, g_v_;
};

struct Errors {
  double u, v;
};

Errors SolveAndMeasure(int n) {
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, n, 0.0, 1.0, n);
  CoupledModel model(mesh);

  mphys::SolverOptions opts;
  opts.tolerance.absolute = 1e-11;
  opts.tolerance.relative = 1e-12;

  mphys::SunContext sunctx;
  mphys::MeshSteadySolver solver(model, opts, sunctx);
  solver.Solve();

  const auto& u = model.fields()[model.u()];
  const auto& v = model.fields()[model.v()];
  Errors e{0.0, 0.0};
  for (int c = 0; c < mesh.NCells(); ++c) {
    const auto& p = mesh.cells[c].centroid;
    e.u = std::max(e.u, std::abs(u[c] - UExact(p[0], p[1])));
    e.v = std::max(e.v, std::abs(v[c] - VExact(p[0], p[1])));
  }
  return e;
}

}  // namespace

// A two-field coupled system solved together must recover the manufactured
// solution and converge at ~2nd order in both fields.
TEST(MeshCoupled, TwoFieldReactionDiffusion2D) {
  const Errors e16 = SolveAndMeasure(16);
  const Errors e32 = SolveAndMeasure(32);

  EXPECT_LT(e32.u, 5e-3) << "u error " << e32.u;
  EXPECT_LT(e32.v, 5e-3) << "v error " << e32.v;

  const double rate_u = std::log(e16.u / e32.u) / std::log(2.0);
  const double rate_v = std::log(e16.v / e32.v) / std::log(2.0);
  EXPECT_GT(rate_u, 1.8) << "u convergence rate " << rate_u;
  EXPECT_GT(rate_v, 1.8) << "v convergence rate " << rate_v;
}
