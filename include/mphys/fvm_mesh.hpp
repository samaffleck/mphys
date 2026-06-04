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
// * Boundary conditions are supplied per patch, indexed by Face::patch:
//     - Dirichlet: the field value on the boundary face.
//     - Neumann:   the outward-normal gradient dphi/dn on the boundary face.
namespace mphys::fvm {

// Per-face gradient in the face reference direction (length NFaces).
//   internal face:  (phi[neighbour] - phi[owner]) / delta
//   boundary face:  derived from the patch boundary condition.
std::vector<double> FaceGrad(const std::vector<double>& phi, const Mesh& mesh,
                             const std::vector<BoundaryCondition>& patch_bcs);

// Conservative divergence of a per-face flux (in the face reference direction)
// → cell-centred result (length NCells):
//   res[owner]     += area * flux / V_owner
//   res[neighbour] -= area * flux / V_neighbour   (internal faces only)
std::vector<double> Div(const std::vector<double>& face_flux, const Mesh& mesh);

// div(D * grad(phi)) with uniform diffusivity.
std::vector<double> Laplacian(const std::vector<double>& phi, double D,
                              const Mesh& mesh,
                              const std::vector<BoundaryCondition>& patch_bcs);

}  // namespace mphys::fvm
