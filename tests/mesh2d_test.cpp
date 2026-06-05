#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/fvm_mesh.hpp"
#include "mphys/topology.hpp"

namespace {

using mphys::DirichletBc;
using mphys::NeumannBc;

std::vector<double> Sample(const mphys::Mesh& mesh, auto fn) {
  std::vector<double> v(mesh.NCells());
  for (int c = 0; c < mesh.NCells(); ++c) {
    const auto& x = mesh.cells[c].centroid;
    v[c] = fn(x[0], x[1]);
  }
  return v;
}

double MaxAbs(const std::vector<double>& v) {
  double m = 0.0;
  for (double x : v) m = std::max(m, std::abs(x));
  return m;
}

}  // namespace

// Structural sanity: cell/face counts, patch sizes, total volume.
TEST(Mesh2D, Structure) {
  constexpr int nx = 5, ny = 3;
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 2.0, nx, 0.0, 1.5, ny);

  EXPECT_EQ(mesh.dim, 2);
  EXPECT_EQ(mesh.NCells(), nx * ny);

  const int internal = (nx - 1) * ny + nx * (ny - 1);
  const int boundary = 2 * ny + 2 * nx;
  EXPECT_EQ(mesh.NFaces(), internal + boundary);

  ASSERT_EQ(mesh.patches.size(), 4u);
  EXPECT_EQ(mesh.patches[0].faces.size(), static_cast<size_t>(ny));  // left
  EXPECT_EQ(mesh.patches[1].faces.size(), static_cast<size_t>(ny));  // right
  EXPECT_EQ(mesh.patches[2].faces.size(), static_cast<size_t>(nx));  // bottom
  EXPECT_EQ(mesh.patches[3].faces.size(), static_cast<size_t>(nx));  // top

  double vol = 0.0;
  for (const auto& c : mesh.cells) vol += c.volume;
  EXPECT_NEAR(vol, 2.0 * 1.5, 1e-12);  // (x1-x0)*(y1-y0)
}

// The face-based Laplacian is exact (to machine precision) for a linear field,
// in both the x- and y-directions, using the unchanged 1D operators.
TEST(Mesh2D, LaplacianExactForLinearField) {
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, 12, 0.0, 1.0, 9);

  // u = 2 + 3x : constant on left/right (Dirichlet), zero y-gradient (Neumann 0).
  {
    auto u = [](double x, double /*y*/) { return 2.0 + 3.0 * x; };
    // patch order: left, right, bottom, top.
    std::vector<mphys::PatchBc> bcs = {
        DirichletBc(2.0 + 3.0 * 0.0), DirichletBc(2.0 + 3.0 * 1.0),
        NeumannBc(0.0), NeumannBc(0.0)};
    const auto lap = mphys::fvm::Laplacian(Sample(mesh, u), 1.0, mesh, bcs);
    EXPECT_LT(MaxAbs(lap), 1e-10);
  }

  // u = 1 - 0.7y : constant on bottom/top (Dirichlet), zero x-gradient (Neumann 0).
  {
    auto u = [](double /*x*/, double y) { return 1.0 - 0.7 * y; };
    std::vector<mphys::PatchBc> bcs = {
        NeumannBc(0.0), NeumannBc(0.0),
        DirichletBc(1.0 - 0.7 * 0.0), DirichletBc(1.0 - 0.7 * 1.0)};
    const auto lap = mphys::fvm::Laplacian(Sample(mesh, u), 1.0, mesh, bcs);
    EXPECT_LT(MaxAbs(lap), 1e-10);
  }
}

// Manufactured solution: u = sin(pi x) sin(pi y) on [0,1]^2 (zero on all edges).
// Analytic Laplacian = -2 pi^2 u. The discrete operator must converge at 2nd
// order in the mesh spacing.
TEST(Mesh2D, LaplacianSecondOrderConvergence) {
  constexpr double pi = std::numbers::pi;
  auto u = [](double x, double y) { return std::sin(pi * x) * std::sin(pi * y); };
  auto lap_exact = [&](double x, double y) { return -2.0 * pi * pi * u(x, y); };

  // All four edges are zero -> Dirichlet 0 everywhere (constant per patch).
  std::vector<mphys::PatchBc> bcs = {
      DirichletBc(0.0), DirichletBc(0.0), DirichletBc(0.0), DirichletBc(0.0)};

  double prev_err = 0.0, prev_h = 0.0;
  double last_rate = 0.0;
  for (int n : {16, 32, 64}) {
    const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, n, 0.0, 1.0, n);
    const auto lap = mphys::fvm::Laplacian(Sample(mesh, u), 1.0, mesh, bcs);

    double err = 0.0;
    for (int c = 0; c < mesh.NCells(); ++c) {
      const auto& xc = mesh.cells[c].centroid;
      err = std::max(err, std::abs(lap[c] - lap_exact(xc[0], xc[1])));
    }

    const double h = 1.0 / n;
    if (prev_err > 0.0) {
      last_rate = std::log(prev_err / err) / std::log(prev_h / h);
    }
    prev_err = err;
    prev_h = h;
  }

  // Expect ~2nd order; allow a little slack for boundary truncation.
  EXPECT_GT(last_rate, 1.9) << "observed convergence rate " << last_rate;
}

// Per-face Dirichlet BCs let a fully 2D-varying boundary be represented exactly.
// A field linear in both x and y has zero Laplacian; with per-face boundary
// values sampled from the exact field, the discrete operator must return ~0 —
// something a single constant value per patch could not achieve.
TEST(Mesh2D, PerFaceDirichletExactForBilinearField) {
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, 10, 0.0, 1.0, 14);
  auto u = [](double x, double y) { return 2.0 + 3.0 * x + 4.0 * y; };

  std::vector<mphys::PatchBc> bcs;
  for (int p = 0; p < 4; ++p) {
    bcs.push_back(mphys::fvm::MakePatchBc(
        mphys::BcType::kDirichlet, mesh, p,
        [&](const mphys::Vec3& c) { return u(c[0], c[1]); }));
  }

  const auto lap = mphys::fvm::Laplacian(Sample(mesh, u), 1.0, mesh, bcs);
  EXPECT_LT(MaxAbs(lap), 1e-10);
}

// A per-face diffusivity that is uniform must reproduce the uniform-D overload.
TEST(Mesh2D, VariableDiffusivityMatchesUniform) {
  constexpr double D = 2.3;
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, 8, 0.0, 1.0, 8);
  auto u = [](double x, double y) { return std::sin(2.0 * x) * (1.0 + y); };
  const auto phi = Sample(mesh, u);

  std::vector<mphys::PatchBc> bcs = {DirichletBc(0.0), DirichletBc(0.0),
                                     DirichletBc(0.0), DirichletBc(0.0)};
  const std::vector<double> D_face(mesh.NFaces(), D);

  const auto uniform = mphys::fvm::Laplacian(phi, D, mesh, bcs);
  const auto variable = mphys::fvm::Laplacian(phi, D_face, mesh, bcs);

  ASSERT_EQ(uniform.size(), variable.size());
  for (size_t i = 0; i < uniform.size(); ++i) {
    EXPECT_NEAR(variable[i], uniform[i], 1e-12) << "cell " << i;
  }
}
