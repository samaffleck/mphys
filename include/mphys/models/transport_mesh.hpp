#pragma once

#include <array>
#include <utility>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/fvm_mesh.hpp"
#include "mphys/mesh_model.hpp"
#include "mphys/topology.hpp"

namespace mphys::models {

// Single-species convection-diffusion-reaction on a face-based Mesh (1D/2D/3D):
//
//   dc/dt + div(u c) - D lap(c) + k c = 0
//
// with a uniform advection velocity u, uniform diffusivity D and a first-order
// reaction rate k. This is the dimension-independent counterpart to the 1D
// ConvDiffReactionModel: the same residual serves both the steady
// (MeshSteadySolver) and transient (MeshTransientSolver) solves — for a steady
// solve `ydot` is empty.
//
// The velocity is given as a Cartesian vector and projected onto each face
// normal once at construction. Boundary conditions are supplied per patch
// (ordered as Mesh::patches) for the single concentration field.
class ConvDiffReactionMesh : public MeshModel {
 public:
  ConvDiffReactionMesh(const Mesh& mesh, std::array<double, 3> velocity,
                       double D, double k, std::vector<PatchBc> bcs)
      : MeshModel(mesh), D_(D), k_(k) {
    c_ = AddField("c", 0.0);
    SetBcs(c_, std::move(bcs));

    // Project the uniform velocity onto each face's reference normal so the
    // upwind Convection operator sees the correct signed face velocity.
    u_face_.resize(mesh.NFaces());
    for (int f = 0; f < mesh.NFaces(); ++f) {
      const auto& n = mesh.faces[f].normal;
      u_face_[f] = velocity[0] * n[0] + velocity[1] * n[1] + velocity[2] * n[2];
    }
  }

  void Residual(double /*t*/, const std::vector<std::vector<double>>& y,
                const std::vector<std::vector<double>>& ydot,
                std::vector<std::vector<double>>& rr) override {
    const auto conv = fvm::Convection(y[c_], u_face_, mesh_, bcs(c_));
    const auto lap = fvm::Laplacian(y[c_], D_, mesh_, bcs(c_));
    const bool transient = !ydot.empty();
    for (int i = 0; i < mesh_.NCells(); ++i) {
      const double rhs = conv[i] - lap[i] + k_ * y[c_][i];
      rr[c_][i] = transient ? ydot[c_][i] + rhs : rhs;
    }
  }

  int c() const { return c_; }

 private:
  int c_ = 0;
  double D_, k_;
  std::vector<double> u_face_;
};

}  // namespace mphys::models
