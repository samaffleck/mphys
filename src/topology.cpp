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
    f.patch_face = 0;
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
    f.patch_face = 0;
    mesh.patches[1].faces.push_back(static_cast<int>(mesh.faces.size()));
    mesh.faces.push_back(f);
  }

  return mesh;
}

Mesh MakeStructuredMesh2D(double x0, double x1, int nx,
                          double y0, double y1, int ny) {
  assert(x1 > x0 && y1 > y0);
  assert(nx > 0 && ny > 0);

  Mesh mesh;
  mesh.dim = 2;
  mesh.coord_system = CoordSystem::kCartesian;

  const double dx = (x1 - x0) / nx;
  const double dy = (y1 - y0) / ny;
  const auto cell_index = [nx](int i, int j) { return j * nx + i; };

  // Cells (row-major: index = j*nx + i).
  mesh.cells.resize(nx * ny);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      Cell& c = mesh.cells[cell_index(i, j)];
      c.centroid = {x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy, 0.0};
      c.volume = dx * dy;
    }
  }

  // Reserve: internal x-faces (nx-1)*ny + internal y-faces nx*(ny-1)
  //          + boundary faces 2*ny + 2*nx.
  mesh.faces.reserve((nx - 1) * ny + nx * (ny - 1) + 2 * ny + 2 * nx);

  // Internal faces normal to x (between (i,j) and (i+1,j)).
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx - 1; ++i) {
      Face f;
      f.owner = cell_index(i, j);
      f.neighbour = cell_index(i + 1, j);
      f.area = dy;
      f.normal = {1.0, 0.0, 0.0};
      f.centroid = {x0 + (i + 1) * dx, y0 + (j + 0.5) * dy, 0.0};
      f.delta = dx;
      mesh.faces.push_back(f);
    }
  }

  // Internal faces normal to y (between (i,j) and (i,j+1)).
  for (int j = 0; j < ny - 1; ++j) {
    for (int i = 0; i < nx; ++i) {
      Face f;
      f.owner = cell_index(i, j);
      f.neighbour = cell_index(i, j + 1);
      f.area = dx;
      f.normal = {0.0, 1.0, 0.0};
      f.centroid = {x0 + (i + 0.5) * dx, y0 + (j + 1) * dy, 0.0};
      f.delta = dy;
      mesh.faces.push_back(f);
    }
  }

  // Boundary patches: left, right, bottom, top.
  mesh.patches.resize(4);
  mesh.patches[0].name = "left";
  mesh.patches[1].name = "right";
  mesh.patches[2].name = "bottom";
  mesh.patches[3].name = "top";

  // Appends one boundary face owned by `owner` to patch `p`.
  const auto add_boundary = [&](int p, int owner, double area, Vec3 normal,
                                Vec3 centroid, double delta) {
    Face f;
    f.owner = owner;
    f.neighbour = -1;
    f.area = area;
    f.normal = normal;
    f.centroid = centroid;
    f.delta = delta;
    f.patch = p;
    f.patch_face = static_cast<int>(mesh.patches[p].faces.size());
    mesh.patches[p].faces.push_back(static_cast<int>(mesh.faces.size()));
    mesh.faces.push_back(f);
  };

  for (int j = 0; j < ny; ++j) {
    const double yc = y0 + (j + 0.5) * dy;
    add_boundary(0, cell_index(0, j), dy, {-1.0, 0.0, 0.0}, {x0, yc, 0.0}, 0.5 * dx);
    add_boundary(1, cell_index(nx - 1, j), dy, {1.0, 0.0, 0.0}, {x1, yc, 0.0}, 0.5 * dx);
  }
  for (int i = 0; i < nx; ++i) {
    const double xc = x0 + (i + 0.5) * dx;
    add_boundary(2, cell_index(i, 0), dx, {0.0, -1.0, 0.0}, {xc, y0, 0.0}, 0.5 * dy);
    add_boundary(3, cell_index(i, ny - 1), dx, {0.0, 1.0, 0.0}, {xc, y1, 0.0}, 0.5 * dy);
  }

  return mesh;
}

Mesh MakeStructuredMesh3D(double x0, double x1, int nx,
                          double y0, double y1, int ny,
                          double z0, double z1, int nz) {
  assert(x1 > x0 && y1 > y0 && z1 > z0);
  assert(nx > 0 && ny > 0 && nz > 0);

  Mesh mesh;
  mesh.dim = 3;
  mesh.coord_system = CoordSystem::kCartesian;

  const double dx = (x1 - x0) / nx;
  const double dy = (y1 - y0) / ny;
  const double dz = (z1 - z0) / nz;
  const auto cell_index = [nx, ny](int i, int j, int k) {
    return (k * ny + j) * nx + i;
  };

  // Cells.
  mesh.cells.resize(nx * ny * nz);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        Cell& c = mesh.cells[cell_index(i, j, k)];
        c.centroid = {x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy,
                      z0 + (k + 0.5) * dz};
        c.volume = dx * dy * dz;
      }
    }
  }

  const double area_x = dy * dz;  // face normal to x
  const double area_y = dx * dz;  // face normal to y
  const double area_z = dx * dy;  // face normal to z

  // Internal faces normal to x.
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx - 1; ++i) {
        Face f;
        f.owner = cell_index(i, j, k);
        f.neighbour = cell_index(i + 1, j, k);
        f.area = area_x;
        f.normal = {1.0, 0.0, 0.0};
        f.centroid = {x0 + (i + 1) * dx, y0 + (j + 0.5) * dy, z0 + (k + 0.5) * dz};
        f.delta = dx;
        mesh.faces.push_back(f);
      }
  // Internal faces normal to y.
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny - 1; ++j)
      for (int i = 0; i < nx; ++i) {
        Face f;
        f.owner = cell_index(i, j, k);
        f.neighbour = cell_index(i, j + 1, k);
        f.area = area_y;
        f.normal = {0.0, 1.0, 0.0};
        f.centroid = {x0 + (i + 0.5) * dx, y0 + (j + 1) * dy, z0 + (k + 0.5) * dz};
        f.delta = dy;
        mesh.faces.push_back(f);
      }
  // Internal faces normal to z.
  for (int k = 0; k < nz - 1; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        Face f;
        f.owner = cell_index(i, j, k);
        f.neighbour = cell_index(i, j, k + 1);
        f.area = area_z;
        f.normal = {0.0, 0.0, 1.0};
        f.centroid = {x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy, z0 + (k + 1) * dz};
        f.delta = dz;
        mesh.faces.push_back(f);
      }

  // Boundary patches: left, right, bottom, top, back, front.
  mesh.patches.resize(6);
  mesh.patches[0].name = "left";
  mesh.patches[1].name = "right";
  mesh.patches[2].name = "bottom";
  mesh.patches[3].name = "top";
  mesh.patches[4].name = "back";
  mesh.patches[5].name = "front";

  const auto add_boundary = [&](int p, int owner, double area, Vec3 normal,
                                Vec3 centroid, double delta) {
    Face f;
    f.owner = owner;
    f.neighbour = -1;
    f.area = area;
    f.normal = normal;
    f.centroid = centroid;
    f.delta = delta;
    f.patch = p;
    f.patch_face = static_cast<int>(mesh.patches[p].faces.size());
    mesh.patches[p].faces.push_back(static_cast<int>(mesh.faces.size()));
    mesh.faces.push_back(f);
  };

  // x-normal boundaries (left x0, right x1).
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j) {
      const double yc = y0 + (j + 0.5) * dy, zc = z0 + (k + 0.5) * dz;
      add_boundary(0, cell_index(0, j, k), area_x, {-1.0, 0.0, 0.0},
                   {x0, yc, zc}, 0.5 * dx);
      add_boundary(1, cell_index(nx - 1, j, k), area_x, {1.0, 0.0, 0.0},
                   {x1, yc, zc}, 0.5 * dx);
    }
  // y-normal boundaries (bottom y0, top y1).
  for (int k = 0; k < nz; ++k)
    for (int i = 0; i < nx; ++i) {
      const double xc = x0 + (i + 0.5) * dx, zc = z0 + (k + 0.5) * dz;
      add_boundary(2, cell_index(i, 0, k), area_y, {0.0, -1.0, 0.0},
                   {xc, y0, zc}, 0.5 * dy);
      add_boundary(3, cell_index(i, ny - 1, k), area_y, {0.0, 1.0, 0.0},
                   {xc, y1, zc}, 0.5 * dy);
    }
  // z-normal boundaries (back z0, front z1).
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const double xc = x0 + (i + 0.5) * dx, yc = y0 + (j + 0.5) * dy;
      add_boundary(4, cell_index(i, j, 0), area_z, {0.0, 0.0, -1.0},
                   {xc, yc, z0}, 0.5 * dz);
      add_boundary(5, cell_index(i, j, nz - 1), area_z, {0.0, 0.0, 1.0},
                   {xc, yc, z1}, 0.5 * dz);
    }

  return mesh;
}

}  // namespace mphys
