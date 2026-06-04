#include "mphys/fvm_mesh.hpp"

#include <cassert>

namespace mphys::fvm {

std::vector<double> FaceGrad(const std::vector<double>& phi, const Mesh& mesh,
                             const std::vector<PatchBc>& patch_bcs) {
  assert(static_cast<int>(phi.size()) == mesh.NCells());

  std::vector<double> g(mesh.NFaces(), 0.0);
  for (int f = 0; f < mesh.NFaces(); ++f) {
    const Face& face = mesh.faces[f];
    if (face.neighbour >= 0) {
      // Internal face: gradient owner -> neighbour.
      g[f] = (phi[face.neighbour] - phi[face.owner]) / face.delta;
    } else {
      // Boundary face: gradient in the outward-normal direction.
      const PatchBc& bc = patch_bcs[face.patch];
      if (bc.type == BcType::kDirichlet) {
        // Ghost-cell Dirichlet: face value at delta from the owner centroid.
        g[f] = (bc.Value(face.patch_face) - phi[face.owner]) / face.delta;
      } else {  // Neumann: value is the outward-normal gradient directly.
        g[f] = bc.Value(face.patch_face);
      }
    }
  }
  return g;
}

std::vector<double> Div(const std::vector<double>& face_flux, const Mesh& mesh) {
  assert(static_cast<int>(face_flux.size()) == mesh.NFaces());

  std::vector<double> res(mesh.NCells(), 0.0);
  for (int f = 0; f < mesh.NFaces(); ++f) {
    const Face& face = mesh.faces[f];
    const double flow = face.area * face_flux[f];
    res[face.owner] += flow / mesh.cells[face.owner].volume;
    if (face.neighbour >= 0) {
      res[face.neighbour] -= flow / mesh.cells[face.neighbour].volume;
    }
  }
  return res;
}

std::vector<double> Laplacian(const std::vector<double>& phi, double D,
                              const Mesh& mesh,
                              const std::vector<PatchBc>& patch_bcs) {
  std::vector<double> flux = FaceGrad(phi, mesh, patch_bcs);
  for (double& v : flux) v *= D;
  return Div(flux, mesh);
}

std::vector<double> Laplacian(const std::vector<double>& phi,
                              const std::vector<double>& D_face, const Mesh& mesh,
                              const std::vector<PatchBc>& patch_bcs) {
  assert(static_cast<int>(D_face.size()) == mesh.NFaces());

  std::vector<double> flux = FaceGrad(phi, mesh, patch_bcs);
  for (int f = 0; f < mesh.NFaces(); ++f) flux[f] *= D_face[f];
  return Div(flux, mesh);
}

std::vector<double> Convection(const std::vector<double>& phi,
                               const std::vector<double>& u_face, const Mesh& mesh,
                               const std::vector<PatchBc>& patch_bcs) {
  assert(static_cast<int>(phi.size()) == mesh.NCells());
  assert(static_cast<int>(u_face.size()) == mesh.NFaces());

  std::vector<double> flux(mesh.NFaces(), 0.0);
  for (int f = 0; f < mesh.NFaces(); ++f) {
    const Face& face = mesh.faces[f];
    const double u = u_face[f];
    double phi_upwind;
    if (face.neighbour >= 0) {
      // Internal face: upwind from the cell the flow comes from.
      phi_upwind = (u >= 0.0) ? phi[face.owner] : phi[face.neighbour];
    } else if (u >= 0.0) {
      // Boundary outflow: transport the interior (owner) value.
      phi_upwind = phi[face.owner];
    } else {
      // Boundary inflow: Dirichlet supplies the incoming value; a Neumann patch
      // is treated as zero-gradient (owner value).
      const PatchBc& bc = patch_bcs[face.patch];
      phi_upwind = (bc.type == BcType::kDirichlet) ? bc.Value(face.patch_face)
                                                   : phi[face.owner];
    }
    flux[f] = u * phi_upwind;
  }
  return Div(flux, mesh);
}

}  // namespace mphys::fvm
