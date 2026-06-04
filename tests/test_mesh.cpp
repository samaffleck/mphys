#include <gtest/gtest.h>

#include "mphys/mesh.hpp"

namespace {

using mphys::CoordSystem;
using mphys::MakeUniformMesh1D;
using mphys::Mesh1D;

TEST(Mesh, ArraySizes) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 5);
  EXPECT_EQ(mesh.n_cells, 5);
  EXPECT_EQ(static_cast<int>(mesh.cell_centres.size()), 5);
  EXPECT_EQ(static_cast<int>(mesh.face_positions.size()), 6);
  EXPECT_EQ(static_cast<int>(mesh.dx.size()), 5);
  EXPECT_EQ(static_cast<int>(mesh.df.size()), 4);  // n_cells - 1
}

TEST(Mesh, UniformSpacing) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 2.0, 4);
  const double dx = 0.5;
  for (int i = 0; i < mesh.n_cells; ++i) {
    EXPECT_NEAR(mesh.dx[i], dx, 1e-15);
    EXPECT_NEAR(mesh.cell_centres[i], (i + 0.5) * dx, 1e-15);
  }
  for (int i = 0; i <= mesh.n_cells; ++i) {
    EXPECT_NEAR(mesh.face_positions[i], i * dx, 1e-15);
  }
}

TEST(Mesh, CentroidDistances) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 4);
  const double dx = 0.25;
  for (double d : mesh.df) EXPECT_NEAR(d, dx, 1e-15);
}

TEST(Mesh, NonZeroOrigin) {
  Mesh1D mesh = MakeUniformMesh1D(1.0, 3.0, 4);
  EXPECT_NEAR(mesh.face_positions.front(), 1.0, 1e-15);
  EXPECT_NEAR(mesh.face_positions.back(), 3.0, 1e-15);
  EXPECT_NEAR(mesh.cell_centres.front(), 1.25, 1e-15);
  EXPECT_NEAR(mesh.cell_centres.back(), 2.75, 1e-15);
}

TEST(Mesh, DefaultCoordSystemIsCartesian) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 3);
  EXPECT_EQ(mesh.coord_system, CoordSystem::kCartesian);
}

TEST(Mesh, ExplicitCoordSystemIsStored) {
  Mesh1D cyl = MakeUniformMesh1D(1.0, 2.0, 3, CoordSystem::kCylindrical);
  Mesh1D sph = MakeUniformMesh1D(1.0, 2.0, 3, CoordSystem::kSpherical);
  EXPECT_EQ(cyl.coord_system, CoordSystem::kCylindrical);
  EXPECT_EQ(sph.coord_system, CoordSystem::kSpherical);
}

TEST(Mesh, SingleCell) {
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 1);
  EXPECT_EQ(mesh.n_cells, 1);
  EXPECT_EQ(static_cast<int>(mesh.cell_centres.size()), 1);
  EXPECT_EQ(static_cast<int>(mesh.face_positions.size()), 2);
  EXPECT_TRUE(mesh.df.empty());
  EXPECT_NEAR(mesh.cell_centres[0], 0.5, 1e-15);
}

}  // namespace
