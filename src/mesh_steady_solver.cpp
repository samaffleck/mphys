#include "mphys/mesh_steady_solver.hpp"

#include <stdexcept>
#include <string>

#include <kinsol/kinsol.h>
#include <kinsol/kinsol_ls.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_iterative.h>

namespace mphys {

MeshSteadySolver::MeshSteadySolver(MeshModel& model, const SolverOptions& opts,
                                   SUNContext sunctx)
    : model_(model),
      opts_(opts),
      n_cells_(model.mesh().NCells()),
      n_fields_(model.NFields()),
      N_(n_cells_ * n_fields_) {
  kin_mem_ = KinMem(sunctx);
  u_ = SunVector(N_, sunctx);

  // Matrix-free GMRES. Use a Krylov subspace large enough to solve small/medium
  // systems without restarts. A Jacobi (diagonal) preconditioner is applied on
  // the right unless disabled.
  const int maxl = std::min(N_, 500);
  const int pretype = opts_.jacobi_preconditioner ? SUN_PREC_RIGHT : SUN_PREC_NONE;
  ls_ = SunLinearSolver(u_, pretype, maxl, sunctx);
  if (opts_.jacobi_preconditioner) precond_.emplace(model_.mesh(), n_fields_);

  // Scratch per-field arrays.
  scratch_y_.assign(n_fields_, std::vector<double>(n_cells_, 0.0));
  scratch_rr_.assign(n_fields_, std::vector<double>(n_cells_, 0.0));

  // Initial guess from the model fields.
  Gather(model_.fields(), u_);

  CheckFlag(KINSetUserData(kin_mem_, this), "KINSetUserData");
  CheckFlag(KINInit(kin_mem_, SystemCb, u_), "KINInit");
  // Matrix-free: pass a NULL matrix so KINSOL uses difference-quotient J*v.
  CheckFlag(KINSetLinearSolver(kin_mem_, ls_, nullptr), "KINSetLinearSolver");
  if (precond_) {
    CheckFlag(KINSetPreconditioner(kin_mem_, PrecSetupCb, PrecSolveCb),
              "KINSetPreconditioner");
  }
  // Constant, moderately tight forcing term: appropriate for the limited
  // accuracy of finite-difference Jacobian-vector products.
  CheckFlag(KINSetEtaForm(kin_mem_, KIN_ETACONSTANT), "KINSetEtaForm");
  CheckFlag(KINSetEtaConstValue(kin_mem_, 1e-4), "KINSetEtaConstValue");
  CheckFlag(KINSetFuncNormTol(kin_mem_, opts_.tolerance.absolute),
            "KINSetFuncNormTol");
  CheckFlag(KINSetScaledStepTol(kin_mem_, opts_.tolerance.relative),
            "KINSetScaledStepTol");
  CheckFlag(KINSetNumMaxIters(kin_mem_, opts_.max_Newton_iter),
            "KINSetNumMaxIters");
  // KINSOL's default maximum Newton step length is derived from ||u0||, which is
  // zero for a zero initial guess and would clamp the first (legitimate) step.
  // Lift the cap so full Newton steps are allowed.
  CheckFlag(KINSetMaxNewtonStep(kin_mem_, 1.0e9), "KINSetMaxNewtonStep");
}

void MeshSteadySolver::Solve() {
  N_Vector u_scale = N_VClone(u_);
  N_Vector f_scale = N_VClone(u_);
  N_VConst(1.0, u_scale);
  N_VConst(1.0, f_scale);

  // Basic Newton (no line search): the systems we target are mild enough that
  // globalisation is unnecessary, and a linear problem converges in one step.
  const int flag = KINSol(kin_mem_, u_, KIN_NONE, u_scale, f_scale);

  N_VDestroy(u_scale);
  N_VDestroy(f_scale);

  if (flag < 0) {
    throw std::runtime_error("KINSol failed with flag " + std::to_string(flag));
  }

  Scatter(u_, model_.fields());
}

// static
int MeshSteadySolver::SystemCb(N_Vector uu, N_Vector fval, void* user_data) {
  auto& solver = *static_cast<MeshSteadySolver*>(user_data);
  solver.Scatter(uu, solver.scratch_y_);
  // Steady solve: no time derivative — pass an empty ydot.
  static const std::vector<std::vector<double>> kNoYdot;
  solver.model_.Residual(0.0, solver.scratch_y_, kNoYdot, solver.scratch_rr_);
  solver.Gather(solver.scratch_rr_, fval);
  return 0;
}

long MeshSteadySolver::NumLinearIterations() const {
  long n = 0;
  KINGetNumLinIters(kin_mem_, &n);
  return n;
}

// static
int MeshSteadySolver::PrecSetupCb(N_Vector u, N_Vector /*uscale*/,
                                  N_Vector fval, N_Vector /*fscale*/,
                                  void* user_data) {
  auto& solver = *static_cast<MeshSteadySolver*>(user_data);
  const double* base_F = N_VGetArrayPointer(fval);

  auto eval = [&solver](const double* y, const double* /*ydot*/, double* F) {
    static const std::vector<std::vector<double>> kNoYdot;
    for (int k = 0; k < solver.n_fields_; ++k) {
      const int base = k * solver.n_cells_;
      for (int i = 0; i < solver.n_cells_; ++i) solver.scratch_y_[k][i] = y[base + i];
    }
    solver.model_.Residual(0.0, solver.scratch_y_, kNoYdot, solver.scratch_rr_);
    for (int k = 0; k < solver.n_fields_; ++k) {
      const int base = k * solver.n_cells_;
      for (int i = 0; i < solver.n_cells_; ++i) F[base + i] = solver.scratch_rr_[k][i];
    }
  };

  solver.precond_->Update(eval, N_VGetArrayPointer(u), nullptr, 0.0, base_F);
  return 0;
}

// static
int MeshSteadySolver::PrecSolveCb(N_Vector /*u*/, N_Vector /*uscale*/,
                                  N_Vector /*fval*/, N_Vector /*fscale*/,
                                  N_Vector v, void* user_data) {
  auto& solver = *static_cast<MeshSteadySolver*>(user_data);
  double* r = N_VGetArrayPointer(v);
  solver.precond_->Apply(r, r);  // in place
  return 0;
}

void MeshSteadySolver::Scatter(N_Vector nv,
                               std::vector<std::vector<double>>& fields) const {
  const double* data = N_VGetArrayPointer(nv);
  for (int k = 0; k < n_fields_; ++k) {
    const int base = k * n_cells_;
    for (int i = 0; i < n_cells_; ++i) fields[k][i] = data[base + i];
  }
}

void MeshSteadySolver::Gather(const std::vector<std::vector<double>>& fields,
                              N_Vector nv) const {
  double* data = N_VGetArrayPointer(nv);
  for (int k = 0; k < n_fields_; ++k) {
    const int base = k * n_cells_;
    for (int i = 0; i < n_cells_; ++i) data[base + i] = fields[k][i];
  }
}

// static
void MeshSteadySolver::CheckFlag(int flag, const char* func_name) {
  if (flag < 0) {
    throw std::runtime_error(std::string(func_name) + " failed with flag " +
                             std::to_string(flag));
  }
}

}  // namespace mphys
