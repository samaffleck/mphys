#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/fvm_mesh.hpp"
#include "mphys/fvm_operators.hpp"
#include "mphys/mesh.hpp"
#include "mphys/topology.hpp"

namespace {

constexpr double kTol = 1e-12;

// Map a legacy FieldBcs (left/right, Neumann expressed as +x gradient) to the
// per-patch boundary conditions the face-based operators expect (Neumann
// expressed as outward-normal gradient). The left outward normal is -x, so a
// left Neumann value flips sign; the right outward normal is +x, so it does not.
std::vector<mphys::PatchBc> ToPatchBcs(const mphys::FieldBcs& bcs) {
  mphys::BoundaryCondition left = bcs.left;
  if (left.type == mphys::BcType::kNeumann) left.value = -left.value;
  return {left, bcs.right};
}

// Sample an analytic profile at the (shared) cell centres of both mesh types.
std::vector<double> Sample(const mphys::Mesh1D& mesh, auto fn) {
  std::vector<double> v(mesh.n_cells);
  for (int i = 0; i < mesh.n_cells; ++i) v[i] = fn(mesh.cell_centres[i]);
  return v;
}

// Assert the face-based Laplacian reproduces the legacy 1D Laplacian cell-by-cell.
void CheckLaplacian(mphys::CoordSystem coord, const mphys::FieldBcs& bcs) {
  constexpr double x0 = 0.5, x1 = 2.5, D = 1.7;
  constexpr int n = 16;

  const mphys::Mesh1D legacy = mphys::MakeUniformMesh1D(x0, x1, n, coord);
  const mphys::Mesh mesh = mphys::MakeMesh1D(x0, x1, n, coord);

  // A non-trivial, non-linear profile so every term is exercised.
  auto profile = [](double x) { return std::sin(1.3 * x) + 0.5 * x * x; };
  const std::vector<double> phi_vals = Sample(legacy, profile);

  mphys::Field phi("phi", n);
  phi.values = phi_vals;

  const mphys::Field legacy_lap = mphys::fvm::Laplacian(phi, D, legacy, bcs);
  const std::vector<double> mesh_lap =
      mphys::fvm::Laplacian(phi_vals, D, mesh, ToPatchBcs(bcs));

  ASSERT_EQ(static_cast<int>(mesh_lap.size()), n);
  for (int i = 0; i < n; ++i) {
    EXPECT_NEAR(mesh_lap[i], legacy_lap[i], kTol) << "coord=" << static_cast<int>(coord)
                                                  << " cell " << i;
  }
}

}  // namespace

TEST(FvmMesh, LaplacianDirichletDirichlet) {
  const mphys::FieldBcs bcs{mphys::DirichletBc(1.0), mphys::DirichletBc(-0.5)};
  CheckLaplacian(mphys::CoordSystem::kCartesian, bcs);
  CheckLaplacian(mphys::CoordSystem::kCylindrical, bcs);
  CheckLaplacian(mphys::CoordSystem::kSpherical, bcs);
}

TEST(FvmMesh, LaplacianDirichletNeumann) {
  // Right boundary Neumann (outward normal == +x, so no sign change).
  const mphys::FieldBcs bcs{mphys::DirichletBc(1.0), mphys::NeumannBc(0.3)};
  CheckLaplacian(mphys::CoordSystem::kCartesian, bcs);
  CheckLaplacian(mphys::CoordSystem::kCylindrical, bcs);
  CheckLaplacian(mphys::CoordSystem::kSpherical, bcs);
}

TEST(FvmMesh, LaplacianNeumannDirichlet) {
  // Left boundary Neumann (outward normal == -x, exercises the sign mapping).
  const mphys::FieldBcs bcs{mphys::NeumannBc(0.4), mphys::DirichletBc(2.0)};
  CheckLaplacian(mphys::CoordSystem::kCartesian, bcs);
  CheckLaplacian(mphys::CoordSystem::kCylindrical, bcs);
  CheckLaplacian(mphys::CoordSystem::kSpherical, bcs);
}

// Face-based upwind convection must reproduce the legacy 1D upwind scheme on
// interior cells (boundary cells use a deliberately different, more standard
// boundary treatment, so they are excluded).
TEST(FvmMesh, ConvectionInteriorMatchesLegacy) {
  constexpr double x0 = 0.0, x1 = 1.0, u = 0.7;
  constexpr int n = 24;

  const mphys::Mesh1D legacy = mphys::MakeUniformMesh1D(x0, x1, n);
  const mphys::Mesh mesh = mphys::MakeMesh1D(x0, x1, n);

  auto profile = [](double x) { return std::sin(2.0 * x) + x; };
  const std::vector<double> phi_vals = Sample(legacy, profile);

  mphys::Field phi("phi", n);
  phi.values = phi_vals;
  const mphys::FieldBcs bcs{mphys::DirichletBc(profile(x0)),
                            mphys::DirichletBc(profile(x1))};
  const mphys::Field legacy_conv = mphys::fvm::Convection(phi, u, legacy, bcs);

  // Project the uniform +x velocity onto each face normal.
  std::vector<double> u_face(mesh.NFaces());
  for (int f = 0; f < mesh.NFaces(); ++f) u_face[f] = u * mesh.faces[f].normal[0];

  const std::vector<double> mesh_conv =
      mphys::fvm::Convection(phi_vals, u_face, mesh, ToPatchBcs(bcs));

  for (int i = 1; i < n - 1; ++i) {
    EXPECT_NEAR(mesh_conv[i], legacy_conv[i], kTol) << "interior cell " << i;
  }
}

// Convecting a uniform field with a uniform velocity (incompressible flow) must
// produce zero divergence everywhere when the inflow Dirichlet value matches.
TEST(FvmMesh, ConvectionUniformFieldConserves) {
  constexpr double C = 3.5, u = 1.2;
  const mphys::Mesh mesh = mphys::MakeMesh1D(0.0, 1.0, 10);

  const std::vector<double> phi(mesh.NCells(), C);
  std::vector<double> u_face(mesh.NFaces());
  for (int f = 0; f < mesh.NFaces(); ++f) u_face[f] = u * mesh.faces[f].normal[0];

  const std::vector<mphys::PatchBc> bcs = {mphys::DirichletBc(C), mphys::DirichletBc(C)};
  const std::vector<double> conv = mphys::fvm::Convection(phi, u_face, mesh, bcs);

  double m = 0.0;
  for (double v : conv) m = std::max(m, std::abs(v));
  EXPECT_LT(m, kTol);
}

// Div of a per-face flux must be conservative: with all face fluxes equal in a
// Cartesian mesh, interior cells (equal in/out areas) see zero divergence.
TEST(FvmMesh, DivConservation) {
  const mphys::Mesh mesh = mphys::MakeMesh1D(0.0, 1.0, 8);
  std::vector<double> flux(mesh.NFaces(), 2.5);
  // Boundary faces use outward-normal reference, so flip their sign to make the
  // physical +x flux uniform across the whole domain.
  for (int f = 0; f < mesh.NFaces(); ++f) {
    if (mesh.faces[f].neighbour < 0 && mesh.faces[f].normal[0] < 0) flux[f] = -2.5;
  }
  const std::vector<double> div = mphys::fvm::Div(flux, mesh);
  for (int i = 0; i < mesh.NCells(); ++i) {
    EXPECT_NEAR(div[i], 0.0, kTol) << "cell " << i;
  }
}
