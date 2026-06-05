#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/fvm_mesh.hpp"
#include "mphys/mesh_model.hpp"
#include "mphys/mesh_transient_solver.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/topology.hpp"

namespace {

constexpr double pi = std::numbers::pi;
constexpr double kD = 1.0;

// 2D transient diffusion: du/dt = D laplacian(u) on [0,1]^2, u = 0 on the
// boundary, u(x,y,0) = sin(pi x) sin(pi y). The lowest mode decays in place:
//   u(x,y,t) = sin(pi x) sin(pi y) exp(-2 D pi^2 t).
class TransientDiffusion2D : public mphys::MeshModel {
 public:
  explicit TransientDiffusion2D(const mphys::Mesh& mesh) : MeshModel(mesh) {
    u_ = AddField("u", 0.0);
    SetBcs(u_, std::vector<mphys::PatchBc>(mesh.patches.size(),
                                           mphys::DirichletBc(0.0)));
    // Initial condition.
    auto& u0 = fields_[u_];
    for (int c = 0; c < mesh.NCells(); ++c) {
      const auto& p = mesh.cells[c].centroid;
      u0[c] = std::sin(pi * p[0]) * std::sin(pi * p[1]);
    }
  }

  void Residual(double /*t*/, const std::vector<std::vector<double>>& y,
                const std::vector<std::vector<double>>& ydot,
                std::vector<std::vector<double>>& rr) override {
    const auto lap = mphys::fvm::Laplacian(y[u_], kD, mesh_, bcs(u_));
    for (int c = 0; c < mesh_.NCells(); ++c) rr[u_][c] = ydot[u_][c] - lap[c];
  }

  int u() const { return u_; }

 private:
  int u_ = 0;
};

double SolveError(int n, double t_end) {
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, n, 0.0, 1.0, n);
  TransientDiffusion2D model(mesh);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-4;
  opts.maximum_time_step = 0.01;
  opts.tolerance.relative = 1e-9;
  opts.tolerance.absolute = 1e-11;

  mphys::SunContext sunctx;
  mphys::MeshTransientSolver solver(model, opts, sunctx);
  const std::string msg = solver.Solve(0.0, t_end);
  EXPECT_TRUE(msg.empty()) << msg;

  const double decay = std::exp(-2.0 * kD * pi * pi * t_end);
  const auto& u = model.fields()[model.u()];
  double err = 0.0;
  for (int c = 0; c < mesh.NCells(); ++c) {
    const auto& p = mesh.cells[c].centroid;
    const double exact = std::sin(pi * p[0]) * std::sin(pi * p[1]) * decay;
    err = std::max(err, std::abs(u[c] - exact));
  }
  return err;
}

}  // namespace

// Transient 2D diffusion must match the analytic decaying mode, and the spatial
// discretisation must converge at ~2nd order (time error is held tight by IDA).
TEST(MeshTransient, Diffusion2DAccuracy) {
  const double err = SolveError(32, 0.05);
  EXPECT_LT(err, 3e-3) << "L_inf error = " << err;
}

TEST(MeshTransient, Diffusion2DSecondOrder) {
  const double e16 = SolveError(16, 0.05);
  const double e32 = SolveError(32, 0.05);
  const double rate = std::log(e16 / e32) / std::log(2.0);
  EXPECT_GT(rate, 1.8) << "rate " << rate << " (e16=" << e16 << ", e32=" << e32
                       << ")";
}
