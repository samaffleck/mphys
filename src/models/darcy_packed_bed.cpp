#include "mphys/models/darcy_packed_bed.hpp"

#include <utility>

#include "mphys/fvm_operators.hpp"

namespace mphys::models {

DarcyPackedBedModel::DarcyPackedBedModel(
    const Mesh1D& mesh, StateVector& sv, std::vector<double> kappa,
    std::vector<double> mu, const std::vector<double>& P_init,
    BoundaryCondition P_lbc, BoundaryCondition P_rbc, BoundaryCondition u_lbc,
    BoundaryCondition u_rbc)
    : Model(mesh, sv), kappa_(std::move(kappa)), mu_(std::move(mu)) {
  P_ = AddField("P", 0.0);
  u_ = AddField("u", 0.0);
  auto& Pf = fields()[P_];
  for (int i = 0; i < mesh_.n_cells; ++i) Pf[i] = P_init[i];
  sv_.MarkFieldAlgebraic("u");  // Darcy has no du/dt — u is algebraic
  SetBcs(P_, {P_lbc, P_rbc});
  SetBcs(u_, {u_lbc, u_rbc});
}

void DarcyPackedBedModel::Residual(double, const std::vector<Field>& y,
                                   const std::vector<Field>& ydot,
                                   const std::vector<double>&,
                                   const std::vector<double>&,
                                   std::vector<Field>& rr,
                                   std::vector<double>&) {
  // ── Darcy's Law (algebraic) ───────────────────────────────────────────────
  //   u + (kappa/mu) * dP/dx = 0
  auto grad_P = fvm::Grad(y[P_], mesh_, bcs_[P_]);
  for (int i = 0; i < mesh_.n_cells; ++i) {
    double dPdx = 0.5 * (grad_P[i] + grad_P[i + 1]);  // cell-centre gradient
    rr[u_][i] = y[u_][i] + (kappa_[i] / mu_[i]) * dPdx;
  }
  // ── Ideal-gas mass balance (rho = PM/RT, T constant → rho ∝ P) ────────────
  //   dP/dt + d(P*u)/dx = 0
  auto u_face = fvm::InterpolateToFaces(y[u_], mesh_, bcs_[u_]);
  rr[P_] = fvm::Ddt(ydot[P_])
         + fvm::Convection(y[P_], u_face, mesh_, bcs_[P_]);
}

}  // namespace mphys::models
