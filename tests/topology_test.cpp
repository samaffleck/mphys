#include <gtest/gtest.h>

#include <cmath>

#include "mphys/mesh.hpp"
#include "mphys/topology.hpp"

namespace {

// Reference coordinate-system metrics, independent of the production code, so
// the test pins down the intended geometry rather than echoing the factory.
double RefArea(double r, mphys::CoordSystem c) {
  switch (c) {
    case mphys::CoordSystem::kCylindrical: return r;
    case mphys::CoordSystem::kSpherical:   return r * r;
    default:                               return 1.0;
  }
}

double RefVolume(double rl, double rr, mphys::CoordSystem c) {
  switch (c) {
    case mphys::CoordSystem::kCylindrical: return (rr * rr - rl * rl) / 2.0;
    case mphys::CoordSystem::kSpherical:
      return (rr * rr * rr - rl * rl * rl) / 3.0;
    default:                               return rr - rl;
  }
}

constexpr double kTol = 1e-12;

// MakeMesh1D must reproduce the geometry of the legacy Mesh1D for every
// coordinate system: same cell count, centroids, volumes, internal-face
// distances, and coordinate-weighted face areas.
void CheckEquivalence(mphys::CoordSystem coord) {
  constexpr double x0 = 0.5, x1 = 2.5;  // offset from 0 so cyl/sph metrics differ per cell
  constexpr int n = 16;

  const mphys::Mesh1D legacy = mphys::MakeUniformMesh1D(x0, x1, n, coord);
  const mphys::Mesh mesh = mphys::MakeMesh1D(x0, x1, n, coord);

  ASSERT_EQ(mesh.NCells(), legacy.n_cells);
  EXPECT_EQ(mesh.dim, 1);
  EXPECT_EQ(mesh.coord_system, coord);

  // Cells: centroid (x) and volume.
  for (int i = 0; i < n; ++i) {
    EXPECT_NEAR(mesh.cells[i].centroid[0], legacy.cell_centres[i], kTol) << "cell " << i;
    const double rl = legacy.face_positions[i];
    const double rr = legacy.face_positions[i + 1];
    EXPECT_NEAR(mesh.cells[i].volume, RefVolume(rl, rr, coord), kTol) << "vol " << i;
  }

  // Faces: n-1 internal + 2 boundary.
  EXPECT_EQ(mesh.NFaces(), n + 1);

  int internal = 0, boundary = 0;
  for (const auto& f : mesh.faces) {
    EXPECT_GE(f.owner, 0);
    EXPECT_LT(f.owner, n);
    if (f.neighbour == -1) {
      ++boundary;
      EXPECT_GE(f.patch, 0);
    } else {
      ++internal;
      EXPECT_EQ(f.patch, -1);
      // Owner/neighbour are adjacent and the face sits between their centroids.
      EXPECT_EQ(f.neighbour, f.owner + 1);
      const double df = legacy.df[f.owner];
      EXPECT_NEAR(f.delta, df, kTol) << "delta owner " << f.owner;
      // Internal face area uses the metric at the shared face position.
      const double x_face = legacy.face_positions[f.owner + 1];
      EXPECT_NEAR(f.area, RefArea(x_face, coord), kTol) << "area owner " << f.owner;
    }
  }
  EXPECT_EQ(internal, n - 1);
  EXPECT_EQ(boundary, 2);

  // Patches: "left" and "right", one face each, outward normals.
  ASSERT_EQ(mesh.patches.size(), 2u);
  EXPECT_EQ(mesh.patches[0].name, "left");
  EXPECT_EQ(mesh.patches[1].name, "right");
  ASSERT_EQ(mesh.patches[0].faces.size(), 1u);
  ASSERT_EQ(mesh.patches[1].faces.size(), 1u);

  const mphys::Face& left = mesh.faces[mesh.patches[0].faces[0]];
  EXPECT_EQ(left.owner, 0);
  EXPECT_NEAR(left.normal[0], -1.0, kTol);
  EXPECT_NEAR(left.area, RefArea(x0, coord), kTol);
  EXPECT_NEAR(left.delta, legacy.cell_centres[0] - x0, kTol);

  const mphys::Face& right = mesh.faces[mesh.patches[1].faces[0]];
  EXPECT_EQ(right.owner, n - 1);
  EXPECT_NEAR(right.normal[0], 1.0, kTol);
  EXPECT_NEAR(right.area, RefArea(x1, coord), kTol);
  EXPECT_NEAR(right.delta, x1 - legacy.cell_centres[n - 1], kTol);
}

}  // namespace

TEST(Topology, MatchesLegacyCartesian) {
  CheckEquivalence(mphys::CoordSystem::kCartesian);
}

TEST(Topology, MatchesLegacyCylindrical) {
  CheckEquivalence(mphys::CoordSystem::kCylindrical);
}

TEST(Topology, MatchesLegacySpherical) {
  CheckEquivalence(mphys::CoordSystem::kSpherical);
}
