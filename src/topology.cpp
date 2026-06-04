#include "mphys/topology.hpp"

#include <cassert>
#include <cmath>

namespace mphys {

namespace {

// Coordinate-system face-area metric at radius/position r.
//   Cartesian: 1,  cylindrical: r,  spherical: r²
double FaceMetric(double r, CoordSystem coord) {
  switch (coord) {
    case CoordSystem::kCylindrical:
      return r;
    case CoordSystem::kSpherical:
      return r * r;
    default:
      return 1.0;
  }
}

// Coordinate-system cell volume between faces at r_l and r_r.
//   Cartesian: r_r-r_l,  cylindrical: (r_r²-r_l²)/2,  spherical: (r_r³-r_l³)/3
double CellVolume(double r_l, double r_r, CoordSystem coord) {
  switch (coord) {
    case CoordSystem::kCylindrical:
      return (r_r * r_r - r_l * r_l) / 2.0;
    case CoordSystem::kSpherical:
      return (r_r * r_r * r_r - r_l * r_l * r_l) / 3.0;
    default:
      return r_r - r_l;
  }
}

}  // namespace

Mesh MakeMesh1D(double x0, double x1, int n_cells, CoordSystem coord) {
  assert(x1 > x0);
  assert(n_cells > 0);

  Mesh mesh;
  mesh.dim = 1;
  mesh.coord_system = coord;

  const double cell_width = (x1 - x0) / n_cells;

  // Cells.
  mesh.cells.resize(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    const double r_l = x0 + i * cell_width;
    const double r_r = r_l + cell_width;
    mesh.cells[i].centroid = {r_l + 0.5 * cell_width, 0.0, 0.0};
    mesh.cells[i].volume = CellVolume(r_l, r_r, coord);
  }

  // Faces: internal faces first (owner i -> neighbour i+1), then the two
  // boundary faces, so internal/boundary blocks are contiguous.
  mesh.faces.reserve(n_cells + 1);
  for (int i = 0; i < n_cells - 1; ++i) {
    const double x_face = x0 + (i + 1) * cell_width;
    Face f;
    f.owner = i;
    f.neighbour = i + 1;
    f.area = FaceMetric(x_face, coord);
    f.normal = {1.0, 0.0, 0.0};  // owner -> neighbour, +x
    f.centroid = {x_face, 0.0, 0.0};
    f.delta = mesh.cells[i + 1].centroid[0] - mesh.cells[i].centroid[0];
    f.patch = -1;
    mesh.faces.push_back(f);
  }

  // Boundary patches.
  mesh.patches.resize(2);
  mesh.patches[0].name = "left";
  mesh.patches[1].name = "right";

  // Left boundary face: owner cell 0, outward normal -x.
  {
    Face f;
    f.owner = 0;
    f.neighbour = -1;
    f.area = FaceMetric(x0, coord);
    f.normal = {-1.0, 0.0, 0.0};
    f.centroid = {x0, 0.0, 0.0};
    f.delta = mesh.cells[0].centroid[0] - x0;
    f.patch = 0;
    mesh.patches[0].faces.push_back(static_cast<int>(mesh.faces.size()));
    mesh.faces.push_back(f);
  }

  // Right boundary face: owner cell n-1, outward normal +x.
  {
    Face f;
    f.owner = n_cells - 1;
    f.neighbour = -1;
    f.area = FaceMetric(x1, coord);
    f.normal = {1.0, 0.0, 0.0};
    f.centroid = {x1, 0.0, 0.0};
    f.delta = x1 - mesh.cells[n_cells - 1].centroid[0];
    f.patch = 1;
    mesh.patches[1].faces.push_back(static_cast<int>(mesh.faces.size()));
    mesh.faces.push_back(f);
  }

  return mesh;
}

}  // namespace mphys
