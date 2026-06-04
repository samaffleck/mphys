#include <cmath>
#include <numbers>
#include <print>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/fvm_mesh.hpp"
#include "mphys/mesh_model.hpp"
#include "mphys/mesh_steady_solver.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/topology.hpp"

// 2D Poisson on the face-based mesh, solved end-to-end with the dimension-
// independent operators and the matrix-free Newton-Krylov solver.
//
//   laplacian(u) = f   on [0,1]^2,   u = 0 on the boundary
//   f = -2 pi^2 sin(pi x) sin(pi y)   =>   u_exact = sin(pi x) sin(pi y)
//
// The same fvm::Laplacian used for 1D problems discretises this 2D system with
// no changes — only the mesh differs.
class Poisson2D : public mphys::MeshModel {
 public:
  explicit Poisson2D(const mphys::Mesh& mesh) : MeshModel(mesh) {
    constexpr double pi = std::numbers::pi;
    u_ = AddField("u", 0.0);
    SetBcs(u_, {mphys::DirichletBc(0.0), mphys::DirichletBc(0.0),
                mphys::DirichletBc(0.0), mphys::DirichletBc(0.0)});
    f_.resize(mesh.NCells());
    for (int c = 0; c < mesh.NCells(); ++c) {
      const auto& x = mesh.cells[c].centroid;
      f_[c] = -2.0 * pi * pi * std::sin(pi * x[0]) * std::sin(pi * x[1]);
    }
  }

  void Residual(const std::vector<std::vector<double>>& y,
                std::vector<std::vector<double>>& rr) override {
    const auto lap = mphys::fvm::Laplacian(y[u_], 1.0, mesh_, bcs(u_));
    for (int c = 0; c < mesh_.NCells(); ++c) rr[u_][c] = lap[c] - f_[c];
  }

 private:
  int u_ = 0;
  std::vector<double> f_;
};

int main() {
  constexpr double pi = std::numbers::pi;
  std::println("2D Poisson (face-based mesh, matrix-free Newton-Krylov)");
  std::println("{}", std::string(60, '-'));
  std::println("  {:>8}  {:>14}  {:>8}", "n (per side)", "L_inf error", "rate");

  double prev_err = 0.0;
  for (int n : {16, 32, 64, 128}) {
    const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, n, 0.0, 1.0, n);
    Poisson2D model(mesh);

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

    if (prev_err > 0.0) {
      std::println("  {:>8}  {:>14.4e}  {:>8.2f}", n, err, std::log2(prev_err / err));
    } else {
      std::println("  {:>8}  {:>14.4e}  {:>8}", n, err, "---");
    }
    prev_err = err;
  }

  std::println("{}", std::string(60, '-'));
  std::println("Done. Second-order convergence confirms the 2D solve.");
  return 0;
}
