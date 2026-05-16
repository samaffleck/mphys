#pragma once

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/mesh.hpp"

namespace mphys::fvm {

// Marker: returns ydot unchanged.
// Usage: rr[0] = fvm::Ddt(ydot[0]) - fvm::Laplacian(...);
Field Ddt(const Field& ydot);

// Face-centred gradient: length n_cells+1.
// Ghost cells are derived from boundary conditions.
Field Grad(const Field& phi, const Mesh1D& mesh, const FieldBcs& bcs);

// Divergence of a face-centred flux field → cell-centred result.
// flux must have length n_cells+1.
Field Div(const Field& face_flux, const Mesh1D& mesh);

// div(D * grad(phi)) with uniform diffusivity.
Field Laplacian(const Field& phi, double D, const Mesh1D& mesh, const FieldBcs& bcs);

// div(D_face * grad(phi)) with spatially varying face diffusivity.
// D_face must have length n_cells+1.
Field Laplacian(const Field& phi, const Field& D_face, const Mesh1D& mesh,
                const FieldBcs& bcs);

// div(u * phi) using first-order upwind scheme with uniform velocity.
Field Convection(const Field& phi, double u, const Mesh1D& mesh, const FieldBcs& bcs);

// div(u_face * phi) using first-order upwind with spatially varying face velocity.
// u_face must have length n_cells+1.
Field Convection(const Field& phi, const Field& u_face, const Mesh1D& mesh,
                 const FieldBcs& bcs);

// Arithmetic-mean interpolation of a cell-centred field to faces.
// Returns a field of length n_cells+1.
Field InterpolateToFaces(const Field& phi, const Mesh1D& mesh, const FieldBcs& bcs);

}  // namespace mphys::fvm
