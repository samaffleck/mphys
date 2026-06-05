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

// Structural sanity for the 3D mesh: counts, patch sizes, total volume.
TEST(Mesh3D, Structure) {
  constexpr int nx = 4, ny = 3, nz = 2;
  const mphys::Mesh mesh =
      mphys::MakeStructuredMesh3D(0.0, 2.0, nx, 0.0, 1.5, ny, 0.0, 1.0, nz);

  EXPECT_EQ(mesh.dim, 3);
  EXPECT_EQ(mesh.NCells(), nx * ny * nz);

  const int internal = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1);
  const int boundary = 2 * (ny * nz) + 2 * (nx * nz) + 2 * (nx * ny);
  EXPECT_EQ(mesh.NFaces(), internal + boundary);

  ASSERT_EQ(mesh.patches.size(), 6u);
  EXPECT_EQ(mesh.patches[0].faces.size(), static_cast<size_t>(ny * nz));  // left
  EXPECT_EQ(mesh.patches[1].faces.size(), static_cast<size_t>(ny * nz));  // right
  EXPECT_EQ(mesh.patches[2].faces.size(), static_cast<size_t>(nx * nz));  // bottom
  EXPECT_EQ(mesh.patches[3].faces.size(), static_cast<size_t>(nx * nz));  // top
  EXPECT_EQ(mesh.patches[4].faces.size(), static_cast<size_t>(nx * ny));  // back
  EXPECT_EQ(mesh.patches[5].faces.size(), static_cast<size_t>(nx * ny));  // front

  double vol = 0.0;
  for (const auto& c : mesh.cells) vol += c.volume;
  EXPECT_NEAR(vol, 2.0 * 1.5 * 1.0, 1e-12);

  // Every boundary face is owned by a valid cell and carries an outward normal.
  for (const auto& f : mesh.faces) {
    if (f.neighbour < 0) {
      EXPECT_GE(f.owner, 0);
      EXPECT_LT(f.owner, mesh.NCells());
      const double nlen = std::sqrt(f.normal[0] * f.normal[0] +
                                    f.normal[1] * f.normal[1] +
                                    f.normal[2] * f.normal[2]);
      EXPECT_NEAR(nlen, 1.0, 1e-12);
    }
  }
}

// 3D Poisson: laplacian(u) = f on [0,1]^3, u = 0 on the boundary.
//   u_exact = sin(pi x) sin(pi y) sin(pi z),  f = -3 pi^2 u_exact.
class Poisson3D : public mphys::MeshModel {
 public:
  explicit Poisson3D(const mphys::Mesh& mesh) : MeshModel(mesh) {
    u_ = AddField("u", 0.0);
    SetBcs(u_, std::vector<mphys::PatchBc>(mesh.patches.size(),
                                           mphys::DirichletBc(0.0)));
    f_.resize(mesh.NCells());
    for (int c = 0; c < mesh.NCells(); ++c) {
      const auto& p = mesh.cells[c].centroid;
      f_[c] = -3.0 * pi * pi * std::sin(pi * p[0]) * std::sin(pi * p[1]) *
              std::sin(pi * p[2]);
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

double SolveError3D(int n) {
  const mphys::Mesh mesh =
      mphys::MakeStructuredMesh3D(0.0, 1.0, n, 0.0, 1.0, n, 0.0, 1.0, n);
  Poisson3D model(mesh);

  mphys::SolverOptions opts;
  opts.tolerance.absolute = 1e-10;
  opts.tolerance.relative = 1e-11;

  mphys::SunContext sunctx;
  mphys::MeshSteadySolver solver(model, opts, sunctx);
  solver.Solve();

  const auto& u = model.fields()[0];
  double err = 0.0;
  for (int c = 0; c < mesh.NCells(); ++c) {
    const auto& p = mesh.cells[c].centroid;
    const double exact =
        std::sin(pi * p[0]) * std::sin(pi * p[1]) * std::sin(pi * p[2]);
    err = std::max(err, std::abs(u[c] - exact));
  }
  return err;
}

}  // namespace

// End-to-end 3D solve with the same operators and solver used in 1D/2D, at
// ~2nd-order accuracy — completing the 1D/2D/3D dimensional story.
TEST(Mesh3D, Poisson3DSecondOrder) {
  const double e8 = SolveError3D(8);
  const double e16 = SolveError3D(16);

  EXPECT_LT(e16, 1e-2) << "L_inf error (16^3) = " << e16;
  const double rate = std::log(e8 / e16) / std::log(2.0);
  EXPECT_GT(rate, 1.8) << "observed convergence rate " << rate
                       << " (e8=" << e8 << ", e16=" << e16 << ")";
}
