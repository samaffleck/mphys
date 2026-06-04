#include "mphys/fvm_mesh.hpp"

#include <cassert>

namespace mphys::fvm {

std::vector<double> FaceGrad(const std::vector<double>& phi, const Mesh& mesh,
                             const std::vector<BoundaryCondition>& patch_bcs) {
  assert(static_cast<int>(phi.size()) == mesh.NCells());

  std::vector<double> g(mesh.NFaces(), 0.0);
  for (int f = 0; f < mesh.NFaces(); ++f) {
    const Face& face = mesh.faces[f];
    if (face.neighbour >= 0) {
      // Internal face: gradient owner -> neighbour.
      g[f] = (phi[face.neighbour] - phi[face.owner]) / face.delta;
    } else {
      // Boundary face: gradient in the outward-normal direction.
      const BoundaryCondition& bc = patch_bcs[face.patch];
      if (bc.type == BcType::kDirichlet) {
        // Ghost-cell Dirichlet: face value is bc.value at delta from the owner.
        g[f] = (bc.value - phi[face.owner]) / face.delta;
      } else {  // Neumann: bc.value is the outward-normal gradient directly.
        g[f] = bc.value;
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
                              const std::vector<BoundaryCondition>& patch_bcs) {
  std::vector<double> flux = FaceGrad(phi, mesh, patch_bcs);
  for (double& v : flux) v *= D;
  return Div(flux, mesh);
}

}  // namespace mphys::fvm
