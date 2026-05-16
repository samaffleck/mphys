#include "mphys/mesh.hpp"

#include <cassert>

namespace mphys {

Mesh1D MakeUniformMesh1D(double x0, double x1, int n_cells, CoordSystem coord) {
  assert(x1 > x0);
  assert(n_cells > 0);

  Mesh1D mesh;
  mesh.n_cells = n_cells;
  mesh.coord_system = coord;

  const double cell_width = (x1 - x0) / n_cells;

  mesh.face_positions.resize(n_cells + 1);
  for (int i = 0; i <= n_cells; ++i) {
    mesh.face_positions[i] = x0 + i * cell_width;
  }

  mesh.cell_centres.resize(n_cells);
  mesh.dx.resize(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    mesh.cell_centres[i] = mesh.face_positions[i] + 0.5 * cell_width;
    mesh.dx[i] = cell_width;
  }

  // df[i] = distance between cell centre i and cell centre i+1 (length n_cells-1)
  mesh.df.resize(n_cells - 1);
  for (int i = 0; i < n_cells - 1; ++i) {
    mesh.df[i] = mesh.cell_centres[i + 1] - mesh.cell_centres[i];
  }

  return mesh;
}

}  // namespace mphys
