#include "mphys/fvm_operators.hpp"

#include <cassert>
#include <cmath>

namespace mphys::fvm {

namespace {

// Ghost cell value and position for the left boundary.
void LeftGhost(const Field& phi, const Mesh1D& mesh, const FieldBcs& bcs,
               double& phi_ghost, double& x_ghost) {
  x_ghost = mesh.cell_centres[0] - mesh.dx[0];
  if (bcs.left.type == BcType::kDirichlet) {
    phi_ghost = 2.0 * bcs.left.value - phi[0];
  } else {  // Neumann: ghost such that (phi[0] - ghost) / dx = flux value
    phi_ghost = phi[0] - bcs.left.value * mesh.dx[0];
  }
}

// Ghost cell value and position for the right boundary.
void RightGhost(const Field& phi, const Mesh1D& mesh, const FieldBcs& bcs,
                double& phi_ghost, double& x_ghost) {
  const int n = phi.NCells();
  x_ghost = mesh.cell_centres[n - 1] + mesh.dx[n - 1];
  if (bcs.right.type == BcType::kDirichlet) {
    phi_ghost = 2.0 * bcs.right.value - phi[n - 1];
  } else {
    phi_ghost = phi[n - 1] + bcs.right.value * mesh.dx[n - 1];
  }
}

// Geometric factors for cell i in the chosen coordinate system.
//
// face area A(r):  1  (Cartesian),  r  (cylindrical),  r²  (spherical)
// cell volume V:   dx, (r_r²-r_l²)/2, (r_r³-r_l³)/3
//
// The FVM flux balance for any coordinate system is:
//   cell residual = (A_right * flux_right - A_left * flux_left) / V
struct GeomFactors {
  double area_left;
  double area_right;
  double volume;
};

static GeomFactors CellGeom(int i, const Mesh1D& mesh) {
  const double r_l = mesh.face_positions[i];
  const double r_r = mesh.face_positions[i + 1];
  switch (mesh.coord_system) {
    case CoordSystem::kCylindrical:
      return {r_l, r_r, (r_r * r_r - r_l * r_l) / 2.0};
    case CoordSystem::kSpherical:
      return {r_l * r_l, r_r * r_r,
              (r_r * r_r * r_r - r_l * r_l * r_l) / 3.0};
    default:  // kCartesian
      return {1.0, 1.0, mesh.dx[i]};
  }
}

}  // namespace

Field Ddt(const Field& ydot) { return ydot; }

Field Grad(const Field& phi, const Mesh1D& mesh, const FieldBcs& bcs) {
  const int n = phi.NCells();
  assert(n == mesh.n_cells);

  Field g("grad_" + phi.name, n + 1, 0.0);

  // Left boundary face
  double phi_ghost, x_ghost;
  LeftGhost(phi, mesh, bcs, phi_ghost, x_ghost);
  g[0] = (phi[0] - phi_ghost) / (mesh.cell_centres[0] - x_ghost);

  // Interior faces
  for (int i = 1; i < n; ++i) {
    g[i] = (phi[i] - phi[i - 1]) / mesh.df[i - 1];
  }

  // Right boundary face
  RightGhost(phi, mesh, bcs, phi_ghost, x_ghost);
  g[n] = (phi_ghost - phi[n - 1]) / (x_ghost - mesh.cell_centres[n - 1]);

  return g;
}

Field Div(const Field& face_flux, const Mesh1D& mesh) {
  const int n = mesh.n_cells;
  assert(face_flux.NCells() == n + 1);

  Field d("div_" + face_flux.name, n, 0.0);
  for (int i = 0; i < n; ++i) {
    const auto g = CellGeom(i, mesh);
    d[i] = (g.area_right * face_flux[i + 1] - g.area_left * face_flux[i]) / g.volume;
  }
  return d;
}

Field Laplacian(const Field& phi, double D, const Mesh1D& mesh, const FieldBcs& bcs) {
  const int n = phi.NCells();
  Field lap("lap_" + phi.name, n, 0.0);

  // Compute face-centred gradients (length n+1) and apply coordinate weighting.
  Field g = Grad(phi, mesh, bcs);
  for (int i = 0; i < n; ++i) {
    const auto geom = CellGeom(i, mesh);
    lap[i] = D * (geom.area_right * g[i + 1] - geom.area_left * g[i]) / geom.volume;
  }

  return lap;
}

Field Laplacian(const Field& phi, const Field& D_face, const Mesh1D& mesh,
                const FieldBcs& bcs) {
  const int n = phi.NCells();
  assert(D_face.NCells() == n + 1);

  Field lap("lap_" + phi.name, n, 0.0);

  Field g = Grad(phi, mesh, bcs);
  for (int i = 0; i < n; ++i) {
    const auto geom = CellGeom(i, mesh);
    lap[i] = (geom.area_right * D_face[i + 1] * g[i + 1] -
              geom.area_left  * D_face[i]     * g[i]) / geom.volume;
  }
  return lap;
}

Field Convection(const Field& phi, double u, const Mesh1D& mesh, const FieldBcs& bcs) {
  const int n = phi.NCells();

  // Build face values using upwind scheme
  Field face("conv_face_" + phi.name, n + 1, 0.0);

  // Left face: upwind from ghost or phi[0]
  double phi_ghost_l, x_ghost_l;
  LeftGhost(phi, mesh, bcs, phi_ghost_l, x_ghost_l);
  face[0] = u >= 0.0 ? u * phi_ghost_l : u * phi[0];

  // Interior faces
  for (int i = 1; i < n; ++i) {
    face[i] = u >= 0.0 ? u * phi[i - 1] : u * phi[i];
  }

  // Right face
  double phi_ghost_r, x_ghost_r;
  RightGhost(phi, mesh, bcs, phi_ghost_r, x_ghost_r);
  face[n] = u >= 0.0 ? u * phi[n - 1] : u * phi_ghost_r;

  return Div(face, mesh);
}

Field Convection(const Field& phi, const Field& u_face, const Mesh1D& mesh,
                 const FieldBcs& bcs) {
  const int n = phi.NCells();
  assert(u_face.NCells() == n + 1);

  Field face("conv_face_" + phi.name, n + 1, 0.0);

  double phi_ghost_l, x_ghost_l;
  LeftGhost(phi, mesh, bcs, phi_ghost_l, x_ghost_l);
  face[0] = u_face[0] >= 0.0 ? u_face[0] * phi_ghost_l : u_face[0] * phi[0];

  for (int i = 1; i < n; ++i) {
    face[i] = u_face[i] >= 0.0 ? u_face[i] * phi[i - 1] : u_face[i] * phi[i];
  }

  double phi_ghost_r, x_ghost_r;
  RightGhost(phi, mesh, bcs, phi_ghost_r, x_ghost_r);
  face[n] = u_face[n] >= 0.0 ? u_face[n] * phi[n - 1] : u_face[n] * phi_ghost_r;

  return Div(face, mesh);
}

Field InterpolateToFaces(const Field& phi, const Mesh1D& mesh, const FieldBcs& bcs) {
  const int n = phi.NCells();
  Field face("face_" + phi.name, n + 1, 0.0);

  // Left face
  double phi_ghost_l, x_ghost_l;
  LeftGhost(phi, mesh, bcs, phi_ghost_l, x_ghost_l);
  face[0] = 0.5 * (phi_ghost_l + phi[0]);

  for (int i = 1; i < n; ++i) {
    face[i] = 0.5 * (phi[i - 1] + phi[i]);
  }

  // Right face
  double phi_ghost_r, x_ghost_r;
  RightGhost(phi, mesh, bcs, phi_ghost_r, x_ghost_r);
  face[n] = 0.5 * (phi[n - 1] + phi_ghost_r);

  return face;
}

}  // namespace mphys::fvm
