#pragma once

#include <array>
#include <string>
#include <vector>

#include "mphys/mesh.hpp"

namespace mphys {

// Minimal 3-component vector. For 1D meshes only the x-component is populated;
// y and z stay zero. Kept deliberately small — geometry, not linear algebra.
using Vec3 = std::array<double, 3>;

// A control volume. The `volume` already bakes in the coordinate-system metric
// (dx for Cartesian, (r_r²-r_l²)/2 cylindrical, (r_r³-r_l³)/3 spherical) so that
// FVM operators are purely topological and never branch on coord_system.
struct Cell {
  Vec3 centroid{0.0, 0.0, 0.0};
  double volume = 0.0;
};

// An internal or boundary face shared by `owner` (always valid) and, for
// internal faces, `neighbour`. The unit `normal` points owner -> neighbour.
// `area` bakes in the coordinate-system metric weight (1, r, or r²) just like
// Cell::volume, so a conservative divergence is simply
//   res[owner]    += area * flux / V_owner
//   res[neighbour]-= area * flux / V_neighbour
// `delta` is the owner→neighbour centroid distance for internal faces, or the
// owner-centroid→face distance for boundary faces (used for gradients/BCs).
struct Face {
  int owner = -1;
  int neighbour = -1;  // -1 for a boundary face
  double area = 0.0;
  Vec3 normal{0.0, 0.0, 0.0};
  Vec3 centroid{0.0, 0.0, 0.0};
  double delta = 0.0;
  int patch = -1;       // index into Mesh::patches, or -1 for an internal face
  int patch_face = -1;  // index of this face within its patch (for per-face BCs)
};

// A named group of boundary faces (e.g. "left", "right"). Boundary conditions
// are attached per patch by the model layer.
struct BoundaryPatch {
  std::string name;
  std::vector<int> faces;  // indices into Mesh::faces
};

// Optional uniform-grid description, populated by the structured factories.
// When `valid`, cells are a tensor-product grid laid out with index
// = (k*ny + j)*nx + i, which is exactly the ordering VTK ImageData expects —
// letting the VTK writer emit proper cells without storing explicit vertices.
struct StructuredInfo {
  bool valid = false;
  int nx = 0, ny = 1, nz = 1;
  double x0 = 0.0, y0 = 0.0, z0 = 0.0;
  double dx = 0.0, dy = 0.0, dz = 0.0;
};

// Dimension-independent, face-based unstructured mesh (OpenFOAM-style).
// Internal faces come first in `faces`, followed by boundary faces, but
// operators rely on Face::neighbour == -1 rather than ordering.
struct Mesh {
  int dim = 1;
  CoordSystem coord_system = CoordSystem::kCartesian;
  std::vector<Cell> cells;
  std::vector<Face> faces;
  std::vector<BoundaryPatch> patches;
  StructuredInfo structured;

  int NCells() const { return static_cast<int>(cells.size()); }
  int NFaces() const { return static_cast<int>(faces.size()); }
};

// Factory: 1D mesh over [x0, x1] with n_cells control volumes, emitted in the
// general face-list form. Produces two boundary patches named "left" and
// "right". Reproduces the geometry of MakeUniformMesh1D(x0, x1, n_cells, coord).
Mesh MakeMesh1D(double x0, double x1, int n_cells,
                CoordSystem coord = CoordSystem::kCartesian);

// Factory: uniform 2D Cartesian mesh over [x0, x1] x [y0, y1] with nx x ny
// cells. Cells are stored row-major (index = j * nx + i). Produces four
// boundary patches named "left" (x0), "right" (x1), "bottom" (y0), "top" (y1).
// The same face-based fvm:: operators that run on a 1D mesh run on this mesh
// unchanged.
Mesh MakeStructuredMesh2D(double x0, double x1, int nx,
                          double y0, double y1, int ny);

// Factory: uniform 3D Cartesian mesh over [x0,x1] x [y0,y1] x [z0,z1] with
// nx x ny x nz cells. Cells are stored with index = (k*ny + j)*nx + i. Produces
// six boundary patches named "left" (x0), "right" (x1), "bottom" (y0),
// "top" (y1), "back" (z0), "front" (z1). The same face-based fvm:: operators
// run on this mesh unchanged.
Mesh MakeStructuredMesh3D(double x0, double x1, int nx,
                          double y0, double y1, int ny,
                          double z0, double z1, int nz);

}  // namespace mphys
