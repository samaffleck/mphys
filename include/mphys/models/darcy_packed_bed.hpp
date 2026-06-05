#pragma once

#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/state_vector.hpp"

namespace mphys::models {

// 1D packed bed with Darcy's law + ideal-gas mass balance.
// Two coupled equations kept explicit:
//   (1) Darcy's Law  (algebraic): u + (kappa/mu)*dP/dx = 0
//   (2) Mass balance (transient): dP/dt + d(P*u)/dx = 0
class DarcyPackedBedModel : public Model {
 public:
  DarcyPackedBedModel(const Mesh1D& mesh, StateVector& sv,
                      std::vector<double> kappa,  // per-cell permeability [m^2]
                      std::vector<double> mu,     // per-cell viscosity    [Pa.s]
                      const std::vector<double>& P_init,  // per-cell initial P
                      BoundaryCondition P_lbc, BoundaryCondition P_rbc,
                      BoundaryCondition u_lbc, BoundaryCondition u_rbc);

  void Residual(double t, const std::vector<Field>& y,
                const std::vector<Field>& ydot, const std::vector<double>& y_alg,
                const std::vector<double>& ydot_alg, std::vector<Field>& rr,
                std::vector<double>& rr_alg) override;

 private:
  int P_ = 0, u_ = 1;
  std::vector<double> kappa_, mu_;
};

}  // namespace mphys::models
