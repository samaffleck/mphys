#include "mphys/transient_solver.hpp"

#include <stdexcept>
#include <string>

#include <ida/ida.h>
#include <ida/ida_ls.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunmatrix/sunmatrix_band.h>

namespace mphys {

TransientSolver::TransientSolver(Model& model, const SolverOptions& opts,
                                 SUNContext sunctx)
    : model_(model), opts_(opts) {
  StateVector& sv = model_.state_vector();
  const int N = sv.TotalSize();
  const int bw = sv.NVars() + sv.NAlgebraics();

  // Allocate RAII SUNDIALS objects
  ida_mem_ = IdaMem(sunctx);
  yy_ = SunVector(N, sunctx);
  yp_ = SunVector(N, sunctx);
  id_ = SunVector(N, sunctx);
  A_  = SunMatrix(static_cast<sunindextype>(N),
                  static_cast<sunindextype>(bw),
                  static_cast<sunindextype>(bw), sunctx);
  ls_ = SunLinearSolver(yy_, A_, sunctx);

  // Populate initial conditions into yy_; zero yp_ (IDACalcIC finds consistent yp)
  sv.Gather(model_.fields(), model_.algebraics(), yy_);
  N_VConst(0.0, yp_);

  // Fill the differential/algebraic id vector
  sv.FillIdVector(id_);

  // Allocate scratch buffers once
  sv.AllocateScratch(scratch_y_, scratch_alg_);
  sv.AllocateScratch(scratch_ydot_, scratch_aldot_);
  sv.AllocateScratch(scratch_rr_, scratch_ralg_);

  // IDA initialisation (order matters)
  CheckFlag(IDASetUserData(ida_mem_, this), "IDASetUserData");
  CheckFlag(IDAInit(ida_mem_, ResidualCb, 0.0, yy_, yp_), "IDAInit");
  CheckFlag(IDASStolerances(ida_mem_, opts_.tolerance.relative,
                            opts_.tolerance.absolute),
            "IDASStolerances");
  CheckFlag(IDASetLinearSolver(ida_mem_, ls_, A_), "IDASetLinearSolver");
  CheckFlag(IDASetId(ida_mem_, id_), "IDASetId");
  CheckFlag(IDASetMaxStep(ida_mem_, opts_.maximum_time_step), "IDASetMaxStep");
  CheckFlag(IDASetMinStep(ida_mem_, opts_.minimum_time_step), "IDASetMinStep");
  CheckFlag(IDASetMaxNumSteps(ida_mem_, 100000), "IDASetMaxNumSteps");
}

void TransientSolver::Solve(
    double t0, double t_end,
    std::function<void(double, const std::vector<Field>&,
                       const std::vector<double>&)>
        output_cb) {
  StateVector& sv = model_.state_vector();

  // Compute consistent initial conditions before the first step
  CheckFlag(IDACalcIC(ida_mem_, IDA_YA_YDP_INIT,
                      t0 + opts_.initial_time_step),
            "IDACalcIC");

  // Fire callback with the initial state
  sv.Scatter(yy_, scratch_y_, scratch_alg_);
  if (output_cb) output_cb(t0, scratch_y_, scratch_alg_);

  // Use IDA_ONE_STEP so the callback fires at every accepted internal step,
  // allowing the caller to filter snapshots at arbitrary intervals.
  CheckFlag(IDASetStopTime(ida_mem_, t_end), "IDASetStopTime");

  double t_ret = t0;
  while (t_ret < t_end) {
    const int flag = IDASolve(ida_mem_, t_end, &t_ret, yy_, yp_, IDA_ONE_STEP);
    if (flag < 0) {
      throw std::runtime_error("IDASolve failed with flag " +
                               std::to_string(flag));
    }
    sv.Scatter(yy_, scratch_y_, scratch_alg_);
    if (output_cb) output_cb(t_ret, scratch_y_, scratch_alg_);
    if (flag == IDA_TSTOP_RETURN || t_ret >= t_end) break;
  }
}

// static
int TransientSolver::ResidualCb(sunrealtype t, N_Vector yy, N_Vector yp,
                                 N_Vector rr, void* user_data) {
  auto& solver = *static_cast<TransientSolver*>(user_data);
  StateVector& sv = solver.model_.state_vector();

  sv.Scatter(yy, solver.scratch_y_, solver.scratch_alg_);
  sv.Scatter(yp, solver.scratch_ydot_, solver.scratch_aldot_);

  solver.model_.Residual(t, solver.scratch_y_, solver.scratch_ydot_,
                         solver.scratch_alg_, solver.scratch_aldot_,
                         solver.scratch_rr_, solver.scratch_ralg_);

  sv.Gather(solver.scratch_rr_, solver.scratch_ralg_, rr);
  return 0;
}

// static
void TransientSolver::CheckFlag(int flag, const char* func_name) {
  if (flag < 0) {
    throw std::runtime_error(std::string(func_name) +
                             " failed with flag " + std::to_string(flag));
  }
}

}  // namespace mphys
