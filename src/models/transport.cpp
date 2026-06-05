#include "mphys/models/transport.hpp"

#include <utility>

#include "mphys/fvm_operators.hpp"

namespace mphys::models {

ConvDiffReactionModel::ConvDiffReactionModel(
    const Mesh1D& mesh, StateVector& sv, Field D_face, Field u_face,
    std::vector<double> k_cell, BoundaryCondition left_bc,
    BoundaryCondition right_bc)
    : Model(mesh, sv),
      D_face_(std::move(D_face)),
      u_face_(std::move(u_face)),
      k_cell_(std::move(k_cell)) {
  c_ = AddField("c", 0.0);
  SetBcs(c_, {left_bc, right_bc});
}

void ConvDiffReactionModel::Residual(double, const std::vector<Field>& y,
                                     const std::vector<Field>& ydot,
                                     const std::vector<double>&,
                                     const std::vector<double>&,
                                     std::vector<Field>& rr,
                                     std::vector<double>&) {
  rr[c_] = fvm::Ddt(ydot[c_])
         + fvm::Convection(y[c_], u_face_, mesh_, bcs_[c_])
         - fvm::Laplacian(y[c_], D_face_, mesh_, bcs_[c_]);
  for (int i = 0; i < mesh_.n_cells; ++i)
    rr[c_][i] += y[c_][i] * k_cell_[i];
}

SteadyDiffusionModel::SteadyDiffusionModel(const Mesh1D& mesh, StateVector& sv,
                                           Field D_face,
                                           BoundaryCondition left_bc,
                                           BoundaryCondition right_bc)
    : Model(mesh, sv), D_face_(std::move(D_face)) {
  double init = 0.5 * (left_bc.value + right_bc.value);
  c_ = AddField("c", init);
  SetBcs(c_, {left_bc, right_bc});
}

void SteadyDiffusionModel::Residual(double, const std::vector<Field>& y,
                                    const std::vector<Field>&,
                                    const std::vector<double>&,
                                    const std::vector<double>&,
                                    std::vector<Field>& rr,
                                    std::vector<double>&) {
  rr[c_] = -fvm::Laplacian(y[c_], D_face_, mesh_, bcs_[c_]);
}

}  // namespace mphys::models
