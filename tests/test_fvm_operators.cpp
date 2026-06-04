#include <gtest/gtest.h>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/fvm_operators.hpp"
#include "mphys/mesh.hpp"

namespace {

using mphys::CoordSystem;
using mphys::DirichletBc;
using mphys::Field;
using mphys::FieldBcs;
using mphys::MakeUniformMesh1D;
using mphys::Mesh1D;
using mphys::NeumannBc;
namespace fvm = mphys::fvm;

// Build a cell-centred field from an analytic function of position.
Field SampleField(const Mesh1D& mesh, const std::string& name, auto fn) {
  Field f(name, mesh.n_cells);
  for (int i = 0; i < mesh.n_cells; ++i) f[i] = fn(mesh.cell_centres[i]);
  return f;
}

// ---------------------------------------------------------------------------
// Ddt
// ---------------------------------------------------------------------------

TEST(Fvm, DdtReturnsInputUnchanged) {
  Field ydot("ydot", 3);
  ydot.values = {1.0, 2.0, 3.0};
  Field out = fvm::Ddt(ydot);
  ASSERT_EQ(out.NCells(), 3);
  for (int i = 0; i < 3; ++i) EXPECT_DOUBLE_EQ(out[i], ydot[i]);
}

// ---------------------------------------------------------------------------
// Grad
// ---------------------------------------------------------------------------

// For a linear field phi = a*x + b with Dirichlet BCs equal to the exact
// boundary values, the face gradient is exactly 'a' everywhere (n+1 faces).
TEST(Fvm, GradOfLinearFieldIsConstantSlope) {
  const double a = 2.0, b = 1.0;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 4);
  Field phi = SampleField(mesh, "phi", [&](double x) { return a * x + b; });
  FieldBcs bcs{DirichletBc(b), DirichletBc(a * 1.0 + b)};

  Field g = fvm::Grad(phi, mesh, bcs);
  ASSERT_EQ(g.NCells(), mesh.n_cells + 1);
  for (int i = 0; i <= mesh.n_cells; ++i) EXPECT_NEAR(g[i], a, 1e-12);
}

// A Neumann BC sets the boundary face gradient directly to its value.
TEST(Fvm, GradNeumannSetsBoundaryFaceGradient) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 4);
  Field phi = SampleField(mesh, "phi", [](double x) { return x * x; });
  FieldBcs bcs{NeumannBc(0.5), NeumannBc(-0.3)};

  Field g = fvm::Grad(phi, mesh, bcs);
  EXPECT_NEAR(g.values.front(), 0.5, 1e-12);
  EXPECT_NEAR(g.values.back(), -0.3, 1e-12);
}

// ---------------------------------------------------------------------------
// Div (Cartesian)
// ---------------------------------------------------------------------------

TEST(Fvm, DivOfConstantFluxIsZeroCartesian) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 5);
  Field flux("flux", mesh.n_cells + 1, 4.2);  // constant on every face
  Field d = fvm::Div(flux, mesh);
  ASSERT_EQ(d.NCells(), mesh.n_cells);
  for (int i = 0; i < mesh.n_cells; ++i) EXPECT_NEAR(d[i], 0.0, 1e-12);
}

TEST(Fvm, DivOfLinearFluxIsConstantCartesian) {
  const double slope = 3.0, intercept = 0.7;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 5);
  Field flux("flux", mesh.n_cells + 1);
  for (int i = 0; i <= mesh.n_cells; ++i)
    flux[i] = slope * mesh.face_positions[i] + intercept;

  Field d = fvm::Div(flux, mesh);
  for (int i = 0; i < mesh.n_cells; ++i) EXPECT_NEAR(d[i], slope, 1e-12);
}

// ---------------------------------------------------------------------------
// Div (curvilinear): conservative form (A_r f_r - A_l f_l) / V
// ---------------------------------------------------------------------------

TEST(Fvm, DivConstantFluxCylindrical) {
  const double F = 2.0;
  Mesh1D mesh = MakeUniformMesh1D(1.0, 2.0, 4, CoordSystem::kCylindrical);
  Field flux("flux", mesh.n_cells + 1, F);
  Field d = fvm::Div(flux, mesh);
  for (int i = 0; i < mesh.n_cells; ++i) {
    const double r_l = mesh.face_positions[i];
    const double r_r = mesh.face_positions[i + 1];
    const double expected = (r_r * F - r_l * F) / ((r_r * r_r - r_l * r_l) / 2.0);
    EXPECT_NEAR(d[i], expected, 1e-12);
  }
}

TEST(Fvm, DivConstantFluxSpherical) {
  const double F = 2.0;
  Mesh1D mesh = MakeUniformMesh1D(1.0, 2.0, 4, CoordSystem::kSpherical);
  Field flux("flux", mesh.n_cells + 1, F);
  Field d = fvm::Div(flux, mesh);
  for (int i = 0; i < mesh.n_cells; ++i) {
    const double r_l = mesh.face_positions[i];
    const double r_r = mesh.face_positions[i + 1];
    const double expected = (r_r * r_r * F - r_l * r_l * F) /
                            ((r_r * r_r * r_r - r_l * r_l * r_l) / 3.0);
    EXPECT_NEAR(d[i], expected, 1e-12);
  }
}

// ---------------------------------------------------------------------------
// Laplacian
// ---------------------------------------------------------------------------

// d²/dx²(linear) = 0 exactly on a uniform Cartesian mesh.
TEST(Fvm, LaplacianOfLinearFieldIsZero) {
  const double a = -1.5, b = 4.0, D = 2.0;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 6);
  Field phi = SampleField(mesh, "phi", [&](double x) { return a * x + b; });
  FieldBcs bcs{DirichletBc(b), DirichletBc(a + b)};

  Field lap = fvm::Laplacian(phi, D, mesh, bcs);
  for (int i = 0; i < mesh.n_cells; ++i) EXPECT_NEAR(lap[i], 0.0, 1e-10);
}

// d²/dx²(x²) = 2, so D*Laplacian = 2D — exact in interior cells.
TEST(Fvm, LaplacianOfQuadraticIsExactInInterior) {
  const double D = 2.0;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 6);
  Field phi = SampleField(mesh, "phi", [](double x) { return x * x; });
  FieldBcs bcs{DirichletBc(0.0), DirichletBc(1.0)};

  Field lap = fvm::Laplacian(phi, D, mesh, bcs);
  for (int i = 1; i < mesh.n_cells - 1; ++i) EXPECT_NEAR(lap[i], 2.0 * D, 1e-9);
}

// The varying-diffusivity overload with a constant D_face must equal the
// uniform-diffusivity overload.
TEST(Fvm, LaplacianVaryingDMatchesUniformWhenConstant) {
  const double D = 1.7;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 5);
  Field phi = SampleField(mesh, "phi", [](double x) { return x * x * x; });
  FieldBcs bcs{DirichletBc(0.0), DirichletBc(1.0)};

  Field d_face("D", mesh.n_cells + 1, D);
  Field lap_uniform = fvm::Laplacian(phi, D, mesh, bcs);
  Field lap_varying = fvm::Laplacian(phi, d_face, mesh, bcs);

  for (int i = 0; i < mesh.n_cells; ++i)
    EXPECT_NEAR(lap_uniform[i], lap_varying[i], 1e-12);
}

// ---------------------------------------------------------------------------
// Convection (first-order upwind)
// ---------------------------------------------------------------------------

// Convection of a uniform field whose Dirichlet BC matches it is zero.
TEST(Fvm, ConvectionOfUniformFieldIsZero) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 5);
  Field phi("phi", mesh.n_cells, 3.0);
  FieldBcs bcs{DirichletBc(3.0), DirichletBc(3.0)};

  Field c = fvm::Convection(phi, 2.0, mesh, bcs);
  for (int i = 0; i < mesh.n_cells; ++i) EXPECT_NEAR(c[i], 0.0, 1e-12);
}

// Hand-computed values pin down upwinding for positive velocity.
// Mesh [0,1], n=2: dx=0.5, faces {0,0.5,1}; phi={1,3}, u=+2,
// left Dirichlet 0 (ghost=-1), right Neumann 0 (ghost=3).
TEST(Fvm, ConvectionUpwindPositiveVelocity) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 2);
  Field phi("phi", 2);
  phi.values = {1.0, 3.0};
  FieldBcs bcs{DirichletBc(0.0), NeumannBc(0.0)};

  Field c = fvm::Convection(phi, 2.0, mesh, bcs);
  ASSERT_EQ(c.NCells(), 2);
  EXPECT_NEAR(c[0], 8.0, 1e-12);
  EXPECT_NEAR(c[1], 8.0, 1e-12);
}

// Same configuration with u=-2 upwinds from the right, giving different values.
TEST(Fvm, ConvectionUpwindNegativeVelocity) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 2);
  Field phi("phi", 2);
  phi.values = {1.0, 3.0};
  FieldBcs bcs{DirichletBc(0.0), NeumannBc(0.0)};

  Field c = fvm::Convection(phi, -2.0, mesh, bcs);
  EXPECT_NEAR(c[0], -8.0, 1e-12);
  EXPECT_NEAR(c[1], 0.0, 1e-12);
}

// The varying-velocity overload with constant u_face matches the uniform one.
TEST(Fvm, ConvectionVaryingVelocityMatchesUniformWhenConstant) {
  const double u = 1.3;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 5);
  Field phi = SampleField(mesh, "phi", [](double x) { return x + 1.0; });
  FieldBcs bcs{DirichletBc(1.0), NeumannBc(0.0)};

  Field u_face("u", mesh.n_cells + 1, u);
  Field c_uniform = fvm::Convection(phi, u, mesh, bcs);
  Field c_varying = fvm::Convection(phi, u_face, mesh, bcs);
  for (int i = 0; i < mesh.n_cells; ++i)
    EXPECT_NEAR(c_uniform[i], c_varying[i], 1e-12);
}

// ---------------------------------------------------------------------------
// InterpolateToFaces
// ---------------------------------------------------------------------------

// Mesh [0,1], n=2; phi={1,3}; left Dirichlet 0 (ghost=-1), right Neumann 0
// (ghost=3). Faces: 0.5*(-1+1)=0, 0.5*(1+3)=2, 0.5*(3+3)=3.
TEST(Fvm, InterpolateToFacesArithmeticMean) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 2);
  Field phi("phi", 2);
  phi.values = {1.0, 3.0};
  FieldBcs bcs{DirichletBc(0.0), NeumannBc(0.0)};

  Field face = fvm::InterpolateToFaces(phi, mesh, bcs);
  ASSERT_EQ(face.NCells(), 3);
  EXPECT_NEAR(face[0], 0.0, 1e-12);
  EXPECT_NEAR(face[1], 2.0, 1e-12);
  EXPECT_NEAR(face[2], 3.0, 1e-12);
}

// Interior face interpolation of a linear field recovers the face position value.
TEST(Fvm, InterpolateToFacesInteriorIsExactForLinear) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 5);
  Field phi = SampleField(mesh, "phi", [](double x) { return 2.0 * x + 1.0; });
  FieldBcs bcs{DirichletBc(1.0), DirichletBc(3.0)};

  Field face = fvm::InterpolateToFaces(phi, mesh, bcs);
  for (int i = 1; i < mesh.n_cells; ++i)
    EXPECT_NEAR(face[i], 2.0 * mesh.face_positions[i] + 1.0, 1e-12);
}

}  // namespace
