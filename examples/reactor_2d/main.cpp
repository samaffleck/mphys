#include <cmath>
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
#include "mphys/vtk_writer.hpp"

// 2D coupled convection-diffusion-reaction in a flow channel — a compact
// multiphysics showcase on the face-based mesh.
//
//   species c:  div(u c) - D_c lap(c) + k c        = 0   (consumed by reaction)
//   energy  T:  div(u T) - D_T lap(T) - beta k c   = 0   (heated by reaction)
//
// Plug flow u = (U, 0). Feed enters cold and pure at the left; the exothermic
// reaction consumes c and releases heat into T. The channel walls are cooled
// (T = 0), so the temperature field is genuinely two-dimensional. The two
// fields are solved together in one MeshModel by the matrix-free solver.
class Reactor2D : public mphys::MeshModel {
 public:
  Reactor2D(const mphys::Mesh& mesh, double U, double Dc, double DT, double k,
            double beta)
      : MeshModel(mesh), Dc_(Dc), DT_(DT), k_(k), beta_(beta) {
    c_ = AddField("c", 0.0);
    T_ = AddField("T", 0.0);

    using mphys::DirichletBc;
    using mphys::NeumannBc;
    // Patch order from MakeStructuredMesh2D: left, right, bottom, top.
    SetBcs(c_, {DirichletBc(1.0), NeumannBc(0.0), NeumannBc(0.0), NeumannBc(0.0)});
    SetBcs(T_, {DirichletBc(0.0), NeumannBc(0.0), DirichletBc(0.0), DirichletBc(0.0)});

    // Project the uniform +x velocity onto each face normal.
    u_face_.resize(mesh.NFaces());
    for (int f = 0; f < mesh.NFaces(); ++f) u_face_[f] = U * mesh.faces[f].normal[0];
  }

  void Residual(double /*t*/, const std::vector<std::vector<double>>& y,
                const std::vector<std::vector<double>>& /*ydot*/,
                std::vector<std::vector<double>>& rr) override {
    const auto conv_c = mphys::fvm::Convection(y[c_], u_face_, mesh_, bcs(c_));
    const auto lap_c = mphys::fvm::Laplacian(y[c_], Dc_, mesh_, bcs(c_));
    const auto conv_T = mphys::fvm::Convection(y[T_], u_face_, mesh_, bcs(T_));
    const auto lap_T = mphys::fvm::Laplacian(y[T_], DT_, mesh_, bcs(T_));

    for (int i = 0; i < mesh_.NCells(); ++i) {
      const double rate = k_ * y[c_][i];
      rr[c_][i] = conv_c[i] - lap_c[i] + rate;            // reaction consumes c
      rr[T_][i] = conv_T[i] - lap_T[i] - beta_ * rate;    // reaction heats T
    }
  }

  int c() const { return c_; }
  int T() const { return T_; }

 private:
  int c_ = 0, T_ = 1;
  double Dc_, DT_, k_, beta_;
  std::vector<double> u_face_;
};

int main() {
  constexpr int nx = 80, ny = 32;
  constexpr double L = 1.0, H = 0.4;
  constexpr double U = 1.0, Dc = 0.01, DT = 0.01, k = 2.0, beta = 1.0;

  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, L, nx, 0.0, H, ny);
  Reactor2D model(mesh, U, Dc, DT, k, beta);

  mphys::SolverOptions opts;
  opts.tolerance.absolute = 1e-10;
  opts.tolerance.relative = 1e-11;

  mphys::SunContext sunctx;
  mphys::MeshSteadySolver solver(model, opts, sunctx);
  solver.Solve();

  const auto& c = model.fields()[model.c()];
  const auto& T = model.fields()[model.T()];
  const auto idx = [nx](int i, int j) { return j * nx + i; };

  // Outlet conversion (mean c over the last column) and peak temperature.
  double outlet_c = 0.0;
  for (int j = 0; j < ny; ++j) outlet_c += c[idx(nx - 1, j)];
  outlet_c /= ny;

  double T_max = 0.0;
  double x_hot = 0.0, y_hot = 0.0;
  for (int cell = 0; cell < mesh.NCells(); ++cell) {
    if (T[cell] > T_max) {
      T_max = T[cell];
      x_hot = mesh.cells[cell].centroid[0];
      y_hot = mesh.cells[cell].centroid[1];
    }
  }

  std::println("2D coupled reactor (convection-diffusion-reaction + heat)");
  std::println("{}", std::string(60, '-'));
  std::println("  mesh            : {} x {}  ({} cells, 2 fields)", nx, ny,
               mesh.NCells());
  std::println("  Peclet (U L/Dc) : {:.0f}", U * L / Dc);
  std::println("  inlet c         : 1.000");
  std::println("  outlet mean c   : {:.4f}", outlet_c);
  std::println("  conversion      : {:.1f}%", 100.0 * (1.0 - outlet_c));
  std::println("  peak T          : {:.4f} at (x={:.3f}, y={:.3f})", T_max, x_hot,
               y_hot);

  // Centreline temperature profile (j = ny/2), downsampled.
  std::println("\n  Centreline T(x) along y = {:.3f}:", mesh.cells[idx(0, ny / 2)].centroid[1]);
  std::println("  {:>8}  {:>10}  {:>10}", "x", "c", "T");
  for (int i = 0; i < nx; i += nx / 10) {
    const int cell = idx(i, ny / 2);
    std::println("  {:>8.3f}  {:>10.4f}  {:>10.4f}", mesh.cells[cell].centroid[0],
                 c[cell], T[cell]);
  }

  mphys::WriteVtk("reactor_2d.vti", model);
  std::println("\n  Wrote reactor_2d.vti (fields c, T) — open in ParaView.");

  std::println("{}", std::string(60, '-'));
  std::println("Done. Two coupled transport fields solved together in 2D.");
  return 0;
}
