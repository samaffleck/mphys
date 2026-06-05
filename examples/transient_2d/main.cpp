#include <cmath>
#include <print>
#include <string>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/fvm_mesh.hpp"
#include "mphys/mesh_model.hpp"
#include "mphys/mesh_transient_solver.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/topology.hpp"

// Transient 2D heat diffusion in a square plate, integrated in time by IDA on
// the face-based mesh. The plate starts cold; the left edge is suddenly held
// hot while the other three edges stay cold. Heat diffuses inward and the
// interior warms toward the steady state.
//
//   dT/dt = D laplacian(T)
class HeatPlate2D : public mphys::MeshModel {
 public:
  HeatPlate2D(const mphys::Mesh& mesh, double D) : MeshModel(mesh), D_(D) {
    T_ = AddField("T", 0.0);  // cold start
    using mphys::DirichletBc;
    // Patch order: left, right, bottom, top.
    SetBcs(T_, {DirichletBc(1.0), DirichletBc(0.0), DirichletBc(0.0),
                DirichletBc(0.0)});
  }

  void Residual(double /*t*/, const std::vector<std::vector<double>>& y,
                const std::vector<std::vector<double>>& ydot,
                std::vector<std::vector<double>>& rr) override {
    const auto lap = mphys::fvm::Laplacian(y[T_], D_, mesh_, bcs(T_));
    for (int c = 0; c < mesh_.NCells(); ++c) rr[T_][c] = ydot[T_][c] - lap[c];
  }

  int T() const { return T_; }

 private:
  int T_ = 0;
  double D_;
};

int main() {
  constexpr int nx = 40, ny = 40;
  constexpr double D = 0.1, t_end = 2.0;

  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, nx, 0.0, 1.0, ny);
  HeatPlate2D model(mesh, D);
  const int centre = (ny / 2) * nx + (nx / 2);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-4;
  opts.maximum_time_step = 0.05;
  opts.tolerance.relative = 1e-7;
  opts.tolerance.absolute = 1e-9;

  std::println("Transient 2D heat plate (IDA on the face-based mesh)");
  std::println("{}", std::string(52, '-'));
  std::println("  mesh {} x {}, D = {}, left edge hot at t=0", nx, ny, D);
  std::println("  {:>8}  {:>14}  {:>14}", "t", "T(centre)", "mean T");

  // Print snapshots at roughly evenly spaced times.
  double next_print = 0.0;
  auto report = [&](double t, const std::vector<std::vector<double>>& fields) {
    if (t < next_print && t < t_end) return;
    next_print += 0.2;
    const auto& T = fields[model.T()];
    double mean = 0.0;
    for (double v : T) mean += v;
    mean /= T.size();
    std::println("  {:>8.3f}  {:>14.5f}  {:>14.5f}", t, T[centre], mean);
  };

  mphys::SunContext sunctx;
  mphys::MeshTransientSolver solver(model, opts, sunctx);
  const std::string msg = solver.Solve(0.0, t_end, report);

  std::println("{}", std::string(52, '-'));
  if (msg.empty()) {
    std::println("Done. Interior warms from 0 toward the steady state.");
  } else {
    std::println("Stopped early: {}", msg);
  }
  return 0;
}
