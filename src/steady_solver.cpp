#include "mphys/steady_solver.hpp"

#include <stdexcept>
#include <string>

#include <kinsol/kinsol.h>
#include <kinsol/kinsol_ls.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunmatrix/sunmatrix_band.h>

namespace mphys {

SteadySolver::SteadySolver(Model& model, const SolverOptions& opts,
                           SUNContext sunctx)
    : model_(model), opts_(opts) {
  StateVector& sv = model_.state_vector();
  const int N = sv.TotalSize();
  const int bw = sv.NVars() + sv.NAlgebraics();

  kin_mem_ = KinMem(sunctx);
  u_  = SunVector(N, sunctx);
  A_  = SunMatrix(static_cast<sunindextype>(N),
                  static_cast<sunindextype>(bw),
                  static_cast<sunindextype>(bw), sunctx);
  ls_ = SunLinearSolver(u_, A_, sunctx);

  // Populate initial guess from model fields
  sv.Gather(model_.fields(), model_.algebraics(), u_);

  // Allocate scratch buffers
  sv.AllocateScratch(scratch_y_, scratch_alg_);
  sv.AllocateScratch(scratch_rr_, scratch_ralg_);

  CheckFlag(KINSetUserData(kin_mem_, this), "KINSetUserData");
  CheckFlag(KINInit(kin_mem_, SystemCb, u_), "KINInit");
  CheckFlag(KINSetLinearSolver(kin_mem_, ls_, A_), "KINSetLinearSolver");
  CheckFlag(KINSetFuncNormTol(kin_mem_, opts_.tolerance.absolute),
            "KINSetFuncNormTol");
  CheckFlag(KINSetScaledStepTol(kin_mem_, opts_.tolerance.relative),
            "KINSetScaledStepTol");
  CheckFlag(KINSetNumMaxIters(kin_mem_, opts_.max_Newton_iter),
            "KINSetNumMaxIters");
}

void SteadySolver::Solve(NewtonStrategy strategy) {
  const int kin_strategy =
      strategy == NewtonStrategy::kLineSearch ? KIN_LINESEARCH : KIN_NONE;

  // Clone scale vectors from u_ so they share the same context and size.
  N_Vector u_sc = N_VClone(u_);
  N_Vector f_sc = N_VClone(u_);
  N_VConst(1.0, u_sc);
  N_VConst(1.0, f_sc);

  const int flag = KINSol(kin_mem_, u_, kin_strategy, u_sc, f_sc);

  N_VDestroy(u_sc);
  N_VDestroy(f_sc);

  if (flag < 0) {
    throw std::runtime_error("KINSol failed with flag " + std::to_string(flag));
  }

  // Copy solution back into model fields
  StateVector& sv = model_.state_vector();
  sv.Scatter(u_, model_.fields(), model_.algebraics());
}

// static
int SteadySolver::SystemCb(N_Vector uu, N_Vector fval, void* user_data) {
  auto& solver = *static_cast<SteadySolver*>(user_data);
  StateVector& sv = solver.model_.state_vector();

  sv.Scatter(uu, solver.scratch_y_, solver.scratch_alg_);

  solver.model_.Residual(0.0, solver.scratch_y_, solver.empty_ydot_,
                         solver.scratch_alg_, solver.empty_aldot_,
                         solver.scratch_rr_, solver.scratch_ralg_);

  sv.Gather(solver.scratch_rr_, solver.scratch_ralg_, fval);
  return 0;
}

// static
void SteadySolver::CheckFlag(int flag, const char* func_name) {
  if (flag < 0) {
    throw std::runtime_error(std::string(func_name) +
                             " failed with flag " + std::to_string(flag));
  }
}

}  // namespace mphys
