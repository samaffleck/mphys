#pragma once

#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/state_vector.hpp"

namespace mphys::models {

// Transient convection-diffusion-reaction with spatially-varying coefficients.
//   dc/dt + div(u c) - div(D grad c) + k c = 0
class ConvDiffReactionModel : public Model {
 public:
  ConvDiffReactionModel(const Mesh1D& mesh, StateVector& sv,
                        Field D_face, Field u_face, std::vector<double> k_cell,
                        BoundaryCondition left_bc, BoundaryCondition right_bc);

  void Residual(double t, const std::vector<Field>& y,
                const std::vector<Field>& ydot, const std::vector<double>& y_alg,
                const std::vector<double>& ydot_alg, std::vector<Field>& rr,
                std::vector<double>& rr_alg) override;

 private:
  int c_ = 0;
  Field D_face_, u_face_;
  std::vector<double> k_cell_;
};

// Steady diffusion with spatially-varying face diffusivity.
//   -div(D grad c) = 0
class SteadyDiffusionModel : public Model {
 public:
  SteadyDiffusionModel(const Mesh1D& mesh, StateVector& sv, Field D_face,
                       BoundaryCondition left_bc, BoundaryCondition right_bc);

  void Residual(double t, const std::vector<Field>& y,
                const std::vector<Field>& ydot, const std::vector<double>& y_alg,
                const std::vector<double>& ydot_alg, std::vector<Field>& rr,
                std::vector<double>& rr_alg) override;

 private:
  int c_ = 0;
  Field D_face_;
};

}  // namespace mphys::models
