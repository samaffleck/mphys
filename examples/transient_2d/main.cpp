#include <cmath>
#include <filesystem>
#include <fstream>
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
#include "mphys/vtk_writer.hpp"

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

  // Dump a VTK time series for ParaView animation.
  const std::filesystem::path out_dir = "transient_2d_out";
  std::filesystem::create_directory(out_dir);
  std::vector<std::pair<double, std::string>> frames;

  // Print snapshots and write a VTK frame at roughly evenly spaced times.
  double next_print = 0.0;
  auto report = [&](double t, const std::vector<std::vector<double>>& fields) {
    if (t < next_print && t < t_end) return;
    next_print += 0.2;
    const auto& T = fields[model.T()];
    double mean = 0.0;
    for (double v : T) mean += v;
    mean /= T.size();
    std::println("  {:>8.3f}  {:>14.5f}  {:>14.5f}", t, T[centre], mean);

    const std::string fname =
        "frame_" + std::to_string(frames.size()) + ".vti";
    mphys::WriteVtkImageData((out_dir / fname).string(), mesh, {"T"}, {&T});
    frames.emplace_back(t, fname);
  };

  mphys::SunContext sunctx;
  mphys::MeshTransientSolver solver(model, opts, sunctx);
  const std::string msg = solver.Solve(0.0, t_end, report);

  // Write a ParaView collection (.pvd) tying the frames to their times.
  {
    std::ofstream pvd(out_dir / "transient_2d.pvd");
    pvd << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"1.0\">\n"
        << "  <Collection>\n";
    for (const auto& [t, file] : frames) {
      pvd << "    <DataSet timestep=\"" << t << "\" file=\"" << file << "\"/>\n";
    }
    pvd << "  </Collection>\n</VTKFile>\n";
  }
  std::println("  Wrote {} VTK frames + transient_2d_out/transient_2d.pvd "
               "(open the .pvd in ParaView to animate).",
               frames.size());

  std::println("{}", std::string(52, '-'));
  if (msg.empty()) {
    std::println("Done. Interior warms from 0 toward the steady state.");
  } else {
    std::println("Stopped early: {}", msg);
  }
  return 0;
}
