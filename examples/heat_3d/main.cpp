#include <cmath>
#include <numbers>
#include <print>
#include <string>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/fvm_mesh.hpp"
#include "mphys/mesh_model.hpp"
#include "mphys/mesh_steady_solver.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/topology.hpp"

// Steady 3D heat conduction in a unit cube, solved with the same face-based
// operators and matrix-free Newton-Krylov solver used for 1D and 2D problems —
// only the mesh changes.
//
//   laplacian(T) = q   on [0,1]^3,   T = 0 on the boundary
//   q = -3 pi^2 sin(pi x) sin(pi y) sin(pi z)  =>  T_exact = sin*sin*sin
class Heat3D : public mphys::MeshModel {
 public:
  explicit Heat3D(const mphys::Mesh& mesh) : MeshModel(mesh) {
    constexpr double pi = std::numbers::pi;
    T_ = AddField("T", 0.0);
    SetBcs(T_, std::vector<mphys::PatchBc>(mesh.patches.size(),
                                           mphys::DirichletBc(0.0)));
    q_.resize(mesh.NCells());
    for (int c = 0; c < mesh.NCells(); ++c) {
      const auto& p = mesh.cells[c].centroid;
      q_[c] = -3.0 * pi * pi * std::sin(pi * p[0]) * std::sin(pi * p[1]) *
              std::sin(pi * p[2]);
    }
  }

  void Residual(double /*t*/, const std::vector<std::vector<double>>& y,
                const std::vector<std::vector<double>>& /*ydot*/,
                std::vector<std::vector<double>>& rr) override {
    const auto lap = mphys::fvm::Laplacian(y[T_], 1.0, mesh_, bcs(T_));
    for (int c = 0; c < mesh_.NCells(); ++c) rr[T_][c] = lap[c] - q_[c];
  }

 private:
  int T_ = 0;
  std::vector<double> q_;
};

int main() {
  constexpr double pi = std::numbers::pi;
  std::println("3D heat conduction (face-based mesh, matrix-free Newton-Krylov)");
  std::println("{}", std::string(62, '-'));
  std::println("  {:>10}  {:>10}  {:>14}  {:>8}", "n per side", "cells",
               "L_inf error", "rate");

  double prev_err = 0.0;
  for (int n : {8, 16, 32}) {
    const mphys::Mesh mesh =
        mphys::MakeStructuredMesh3D(0.0, 1.0, n, 0.0, 1.0, n, 0.0, 1.0, n);
    Heat3D model(mesh);

    mphys::SolverOptions opts;
    opts.tolerance.absolute = 1e-10;
    opts.tolerance.relative = 1e-11;

    mphys::SunContext sunctx;
    mphys::MeshSteadySolver solver(model, opts, sunctx);
    solver.Solve();

    const auto& T = model.fields()[0];
    double err = 0.0;
    for (int c = 0; c < mesh.NCells(); ++c) {
      const auto& p = mesh.cells[c].centroid;
      const double exact =
          std::sin(pi * p[0]) * std::sin(pi * p[1]) * std::sin(pi * p[2]);
      err = std::max(err, std::abs(T[c] - exact));
    }

    if (prev_err > 0.0) {
      std::println("  {:>10}  {:>10}  {:>14.4e}  {:>8.2f}", n, mesh.NCells(), err,
                   std::log2(prev_err / err));
    } else {
      std::println("  {:>10}  {:>10}  {:>14.4e}  {:>8}", n, mesh.NCells(), err,
                   "---");
    }
    prev_err = err;
  }

  std::println("{}", std::string(62, '-'));
  std::println("Done. Same operators solve 1D, 2D and 3D — only the mesh differs.");
  return 0;
}
