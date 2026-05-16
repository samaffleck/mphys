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
    d[i] = (face_flux[i + 1] - face_flux[i]) / mesh.dx[i];
  }
  return d;
}

Field Laplacian(const Field& phi, double D, const Mesh1D& mesh, const FieldBcs& bcs) {
  const int n = phi.NCells();
  Field lap("lap_" + phi.name, n, 0.0);

  // Left boundary
  double phi_ghost_l, x_ghost_l;
  LeftGhost(phi, mesh, bcs, phi_ghost_l, x_ghost_l);
  const double grad_left = (phi[0] - phi_ghost_l) / (mesh.cell_centres[0] - x_ghost_l);

  // Right boundary
  double phi_ghost_r, x_ghost_r;
  RightGhost(phi, mesh, bcs, phi_ghost_r, x_ghost_r);
  const double grad_right = (phi_ghost_r - phi[n - 1]) / (x_ghost_r - mesh.cell_centres[n - 1]);

  // Left cell: uses left-boundary face and first interior face
  if (n == 1) {
    lap[0] = D * (grad_right - grad_left) / mesh.dx[0];
  } else {
    const double grad_r0 = (phi[1] - phi[0]) / mesh.df[0];
    lap[0] = D * (grad_r0 - grad_left) / mesh.dx[0];

    // Interior cells
    for (int i = 1; i < n - 1; ++i) {
      const double gr = (phi[i + 1] - phi[i]) / mesh.df[i];
      const double gl = (phi[i] - phi[i - 1]) / mesh.df[i - 1];
      lap[i] = D * (gr - gl) / mesh.dx[i];
    }

    // Right cell
    const double grad_l_last = (phi[n - 1] - phi[n - 2]) / mesh.df[n - 2];
    lap[n - 1] = D * (grad_right - grad_l_last) / mesh.dx[n - 1];
  }

  return lap;
}

Field Laplacian(const Field& phi, const Field& D_face, const Mesh1D& mesh,
                const FieldBcs& bcs) {
  const int n = phi.NCells();
  assert(D_face.NCells() == n + 1);

  Field lap("lap_" + phi.name, n, 0.0);

  // Build all face gradients
  Field g = Grad(phi, mesh, bcs);

  for (int i = 0; i < n; ++i) {
    lap[i] = (D_face[i + 1] * g[i + 1] - D_face[i] * g[i]) / mesh.dx[i];
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
