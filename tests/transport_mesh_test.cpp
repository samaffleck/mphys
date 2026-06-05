#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/mesh_steady_solver.hpp"
#include "mphys/models/transport_mesh.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/topology.hpp"

namespace {

using mphys::DirichletBc;
using mphys::NeumannBc;
using mphys::models::ConvDiffReactionMesh;

mphys::SolverOptions TightOpts() {
  mphys::SolverOptions opts;
  opts.tolerance.absolute = 1e-11;
  opts.tolerance.relative = 1e-12;
  return opts;
}

// ---------------------------------------------------------------------------
// Diffusion-reaction limit (u = 0): D c'' - k c = 0 on [0,1] with c(0)=1,
// c(1)=0 and zero-Neumann on the transverse boundaries. The solution is
// 1D, c(x) = sinh(m(1-x))/sinh(m) with m = sqrt(k/D), so a 2D/3D solve must
// (a) match it and (b) be invariant across the transverse directions.
// ---------------------------------------------------------------------------
double DiffReactExact(double x, double D, double k) {
  const double m = std::sqrt(k / D);
  return std::sinh(m * (1.0 - x)) / std::sinh(m);
}

double DiffReaction2DError(int nx, int ny) {
  constexpr double D = 0.1, k = 1.0;
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, nx, 0.0, 1.0, ny);
  // Patch order: left, right, bottom, top.
  ConvDiffReactionMesh model(
      mesh, {0.0, 0.0, 0.0}, D, k,
      {DirichletBc(1.0), DirichletBc(0.0), NeumannBc(0.0), NeumannBc(0.0)});

  mphys::SunContext sunctx;
  mphys::MeshSteadySolver solver(model, TightOpts(), sunctx);
  solver.Solve();

  const auto& c = model.fields()[model.c()];
  double err = 0.0;
  for (int cell = 0; cell < mesh.NCells(); ++cell) {
    const double exact = DiffReactExact(mesh.cells[cell].centroid[0], D, k);
    err = std::max(err, std::abs(c[cell] - exact));
  }
  return err;
}

TEST(ConvDiffReactionMesh, DiffusionReaction2DMatchesAnalytic) {
  EXPECT_LT(DiffReaction2DError(64, 16), 2e-3);
}

TEST(ConvDiffReactionMesh, DiffusionReaction2DSecondOrder) {
  const double e1 = DiffReaction2DError(16, 16);
  const double e2 = DiffReaction2DError(32, 32);
  const double rate = std::log(e1 / e2) / std::log(2.0);
  EXPECT_GT(rate, 1.8) << "rate " << rate << " (e16=" << e1 << ", e32=" << e2 << ")";
}

// The 2D solution must not vary across y (the equation and BCs are y-uniform).
TEST(ConvDiffReactionMesh, DiffusionReaction2DTransverseInvariance) {
  constexpr int nx = 32, ny = 16;
  constexpr double D = 0.1, k = 1.0;
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, nx, 0.0, 1.0, ny);
  ConvDiffReactionMesh model(
      mesh, {0.0, 0.0, 0.0}, D, k,
      {DirichletBc(1.0), DirichletBc(0.0), NeumannBc(0.0), NeumannBc(0.0)});
  mphys::SunContext sunctx;
  mphys::MeshSteadySolver solver(model, TightOpts(), sunctx);
  solver.Solve();

  const auto& c = model.fields()[model.c()];
  const auto idx = [nx](int i, int j) { return j * nx + i; };
  double spread = 0.0;
  for (int i = 0; i < nx; ++i)
    for (int j = 1; j < ny; ++j)
      spread = std::max(spread, std::abs(c[idx(i, j)] - c[idx(i, 0)]));
  EXPECT_LT(spread, 1e-9) << "transverse spread " << spread;
}

// Same physics in 3D: match the 1D analytic solution and stay invariant across
// both transverse directions.
TEST(ConvDiffReactionMesh, DiffusionReaction3DMatchesAnalytic) {
  constexpr int nx = 32, ny = 8, nz = 8;
  constexpr double D = 0.1, k = 1.0;
  const mphys::Mesh mesh =
      mphys::MakeStructuredMesh3D(0.0, 1.0, nx, 0.0, 1.0, ny, 0.0, 1.0, nz);
  // Patch order: left, right, bottom, top, back, front.
  ConvDiffReactionMesh model(
      mesh, {0.0, 0.0, 0.0}, D, k,
      {DirichletBc(1.0), DirichletBc(0.0), NeumannBc(0.0), NeumannBc(0.0),
       NeumannBc(0.0), NeumannBc(0.0)});
  mphys::SunContext sunctx;
  mphys::MeshSteadySolver solver(model, TightOpts(), sunctx);
  solver.Solve();

  const auto& c = model.fields()[model.c()];
  double err = 0.0;
  for (int cell = 0; cell < mesh.NCells(); ++cell) {
    const double exact = DiffReactExact(mesh.cells[cell].centroid[0], D, k);
    err = std::max(err, std::abs(c[cell] - exact));
  }
  EXPECT_LT(err, 5e-3) << "L_inf error " << err;
}

// ---------------------------------------------------------------------------
// Convection-reaction limit (D -> 0): U c' + k c = 0 with c(0)=1 gives the
// exponential decay c(x) = exp(-(k/U) x). First-order upwind adds numerical
// diffusion, so compare on a fine mesh away from the outflow boundary. This
// pins down the sign/assembly of the Convection operator in 2D.
// ---------------------------------------------------------------------------
TEST(ConvDiffReactionMesh, ConvectionReactionDecaysExponentially) {
  constexpr int nx = 200, ny = 4;
  constexpr double U = 1.0, D = 1e-3, k = 1.0;
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 1.0, nx, 0.0, 1.0, ny);
  // Inflow Dirichlet at left, zero-gradient outflow at right.
  ConvDiffReactionMesh model(
      mesh, {U, 0.0, 0.0}, D, k,
      {DirichletBc(1.0), NeumannBc(0.0), NeumannBc(0.0), NeumannBc(0.0)});
  mphys::SunContext sunctx;
  mphys::MeshSteadySolver solver(model, TightOpts(), sunctx);
  solver.Solve();

  const auto& c = model.fields()[model.c()];
  double err = 0.0;
  for (int cell = 0; cell < mesh.NCells(); ++cell) {
    const double x = mesh.cells[cell].centroid[0];
    if (x > 0.8) continue;  // skip the outflow boundary layer
    err = std::max(err, std::abs(c[cell] - std::exp(-(k / U) * x)));
  }
  EXPECT_LT(err, 0.03) << "L_inf error " << err;
}

}  // namespace
