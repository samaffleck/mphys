#pragma once

#include <vector>

namespace mphys {

enum class CoordSystem { kCartesian, kCylindrical, kSpherical };

struct Mesh1D {
  int n_cells = 0;
  std::vector<double> cell_centres;   // length n_cells
  std::vector<double> face_positions; // length n_cells + 1
  std::vector<double> dx;             // cell widths
  std::vector<double> df;             // centroid-to-centroid distances (for grad reconstruction)
  CoordSystem coord_system = CoordSystem::kCartesian;
};

// Factory: uniform Cartesian mesh from x0 to x1 with n_cells control volumes.
Mesh1D MakeUniformMesh1D(double x0, double x1, int n_cells,
                         CoordSystem coord = CoordSystem::kCartesian);

}  // namespace mphys
