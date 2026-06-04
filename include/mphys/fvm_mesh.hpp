#pragma once

#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/topology.hpp"

// Face-based finite-volume operators over the dimension-independent Mesh.
//
// These are the forward-looking replacements for the 1D-stencil operators in
// fvm_operators.hpp. They loop over faces and scatter contributions to the
// owning cells, so the identical code runs for 1D, 2D and 3D meshes — the only
// difference is the mesh. Coordinate-system metrics (Cartesian/cylindrical/
// spherical) are baked into Face::area and Cell::volume, so nothing here
// branches on the coordinate system.
//
// Conventions
// -----------
// * Cell arrays have length Mesh::NCells(); face arrays have length
//   Mesh::NFaces().
// * A face quantity is expressed in the face's reference direction:
//   owner -> neighbour for an internal face, outward-normal for a boundary face.
// * Boundary conditions are supplied per patch (std::vector<PatchBc> indexed by
//   Face::patch); the per-face value is selected by Face::patch_face, so a patch
//   may carry spatially-varying boundary data.
//     - Dirichlet: the field value on the boundary face.
//     - Neumann:   the outward-normal gradient dphi/dn on the boundary face.
namespace mphys::fvm {

// Per-face gradient in the face reference direction (length NFaces).
//   internal face:  (phi[neighbour] - phi[owner]) / delta
//   boundary face:  derived from the patch boundary condition.
std::vector<double> FaceGrad(const std::vector<double>& phi, const Mesh& mesh,
                             const std::vector<PatchBc>& patch_bcs);

// Conservative divergence of a per-face flux (in the face reference direction)
// → cell-centred result (length NCells):
//   res[owner]     += area * flux / V_owner
//   res[neighbour] -= area * flux / V_neighbour   (internal faces only)
std::vector<double> Div(const std::vector<double>& face_flux, const Mesh& mesh);

// div(D * grad(phi)) with uniform diffusivity.
std::vector<double> Laplacian(const std::vector<double>& phi, double D,
                              const Mesh& mesh,
                              const std::vector<PatchBc>& patch_bcs);

// div(D_face * grad(phi)) with a per-face diffusivity (length NFaces).
std::vector<double> Laplacian(const std::vector<double>& phi,
                              const std::vector<double>& D_face, const Mesh& mesh,
                              const std::vector<PatchBc>& patch_bcs);

// div(u * phi) with first-order upwinding. `u_face` is the velocity component
// along each face's reference direction (length NFaces): owner -> neighbour for
// internal faces, outward-normal for boundary faces. At a boundary the upwind
// value is the owner cell on outflow and, on inflow, the Dirichlet boundary
// value (or the owner cell, i.e. zero-gradient, for a Neumann patch).
std::vector<double> Convection(const std::vector<double>& phi,
                               const std::vector<double>& u_face, const Mesh& mesh,
                               const std::vector<PatchBc>& patch_bcs);

// Build a per-face PatchBc for `patch` by sampling `fn` (a callable taking a
// face centroid Vec3 and returning a double) at every face on the patch, in
// patch-face order. Lets boundary data vary along a 2D/3D boundary.
template <typename Fn>
PatchBc MakePatchBc(BcType type, const Mesh& mesh, int patch, Fn fn) {
  PatchBc bc;
  bc.type = type;
  const auto& faces = mesh.patches[patch].faces;
  bc.values.resize(faces.size());
  for (std::size_t k = 0; k < faces.size(); ++k) {
    bc.values[k] = fn(mesh.faces[faces[k]].centroid);
  }
  return bc;
}

}  // namespace mphys::fvm
