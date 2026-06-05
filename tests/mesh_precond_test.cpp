#include <gtest/gtest.h>

#include <cmath>
#include <utility>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/fvm_mesh.hpp"
#include "mphys/mesh_model.hpp"
#include "mphys/mesh_steady_solver.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/topology.hpp"

namespace {

// Strongly variable-coefficient diffusion: D(x,y) = exp(3x + 3y) spans ~400x
// across the domain, so the Jacobian diagonal varies a lot and a Jacobi
// preconditioner should meaningfully cut GMRES iterations. Driven purely by
// Dirichlet boundaries (no source).
class VarCoeffDiffusion : public mphys::MeshModel {
 public:
  explicit VarCoeffDiffusion(const mphys::Mesh& mesh) : MeshModel(mesh) {
    u_ = AddField("u", 0.5);
    using mphys::DirichletBc;
    SetBcs(u_, {DirichletBc(1.0), DirichletBc(0.0), DirichletBc(0.0),
                DirichletBc(1.0)});
    d_face_.resize(mesh.NFaces());
    for (int f = 0; f < mesh.NFaces(); ++f) {
      const auto& p = mesh.faces[f].centroid;
      d_face_[f] = std::exp(3.0 * p[0] + 3.0 * p[1]);
    }
  }

  void Residual(double /*t*/, const std::vector<std::vector<double>>& y,
                const std::vector<std::vector<double>>& /*ydot*/,
                std::vector<std::vector<double>>& rr) override {
    rr[u_] = mphys::fvm::Laplacian(y[u_], d_face_, mesh_, bcs(u_));
  }

  int u() const { return u_; }

 private:
  int u_ = 0;
  std::vector<double> d_face_;
};

std::pair<long, std::vector<double>> Solve(int n, bool preconditioner) {
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, n, 0.0, 1.0, n);
  VarCoeffDiffusion model(mesh);

  mphys::SolverOptions opts;
  opts.tolerance.absolute = 1e-9;
  opts.tolerance.relative = 1e-10;
  opts.jacobi_preconditioner = preconditioner;

  mphys::SunContext sunctx;
  mphys::MeshSteadySolver solver(model, opts, sunctx);
  solver.Solve();
  return {solver.NumLinearIterations(), model.fields()[model.u()]};
}

}  // namespace

// The Jacobi preconditioner must (a) yield the same solution as the
// unpreconditioned solve and (b) take strictly fewer GMRES iterations on this
// ill-scaled variable-coefficient problem.
TEST(MeshPrecond, ReducesIterationsVariableCoefficient) {
  const auto [iters_off, u_off] = Solve(48, false);
  const auto [iters_on, u_on] = Solve(48, true);

  ASSERT_EQ(u_off.size(), u_on.size());
  double diff = 0.0;
  for (std::size_t i = 0; i < u_off.size(); ++i) {
    diff = std::max(diff, std::abs(u_off[i] - u_on[i]));
  }
  EXPECT_LT(diff, 1e-6) << "preconditioned and unpreconditioned solutions differ";

  EXPECT_GT(iters_off, 0);
  EXPECT_LT(iters_on, iters_off)
      << "linear iterations: preconditioned=" << iters_on
      << " unpreconditioned=" << iters_off;

  std::printf(
      "[ precond ] GMRES iterations: unpreconditioned=%ld, Jacobi=%ld (%.1fx)\n",
      iters_off, iters_on,
      iters_on > 0 ? static_cast<double>(iters_off) / iters_on : 0.0);
}
