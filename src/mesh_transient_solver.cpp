#include "mphys/mesh_transient_solver.hpp"

#include <algorithm>
#include <cstdio>
#include <stdexcept>
#include <string>

#include <ida/ida.h>
#include <ida/ida_ls.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_iterative.h>

namespace mphys {

MeshTransientSolver::MeshTransientSolver(MeshModel& model,
                                         const SolverOptions& opts,
                                         SUNContext sunctx)
    : model_(model),
      opts_(opts),
      n_cells_(model.mesh().NCells()),
      n_fields_(model.NFields()),
      N_(n_cells_ * n_fields_) {
  ida_mem_ = IdaMem(sunctx);
  yy_ = SunVector(N_, sunctx);
  yp_ = SunVector(N_, sunctx);
  id_ = SunVector(N_, sunctx);

  // Matrix-free GMRES, no preconditioner (same approach as the steady solver).
  const int maxl = std::min(N_, 500);
  ls_ = SunLinearSolver(yy_, SUN_PREC_NONE, maxl, sunctx);

  scratch_y_.assign(n_fields_, std::vector<double>(n_cells_, 0.0));
  scratch_ydot_.assign(n_fields_, std::vector<double>(n_cells_, 0.0));
  scratch_rr_.assign(n_fields_, std::vector<double>(n_cells_, 0.0));

  // Initial state from the model fields; start ydot at 0 (IDACalcIC refines it).
  Gather(model_.fields(), yy_);
  N_VConst(0.0, yp_);

  // id vector: 1 for differential components, 0 for algebraic.
  {
    double* id = N_VGetArrayPointer(id_);
    for (int k = 0; k < n_fields_; ++k) {
      const double flag = model_.field_is_differential(k) ? 1.0 : 0.0;
      const int base = k * n_cells_;
      for (int i = 0; i < n_cells_; ++i) id[base + i] = flag;
    }
  }

  CheckFlag(IDASetUserData(ida_mem_, this), "IDASetUserData");
  CheckFlag(IDAInit(ida_mem_, ResidualCb, 0.0, yy_, yp_), "IDAInit");
  CheckFlag(IDASStolerances(ida_mem_, opts_.tolerance.relative,
                            opts_.tolerance.absolute),
            "IDASStolerances");
  CheckFlag(IDASetLinearSolver(ida_mem_, ls_, nullptr), "IDASetLinearSolver");
  CheckFlag(IDASetId(ida_mem_, id_), "IDASetId");
  CheckFlag(IDASetMaxStep(ida_mem_, opts_.maximum_time_step), "IDASetMaxStep");
  CheckFlag(IDASetMaxNumSteps(ida_mem_, 100000), "IDASetMaxNumSteps");
}

std::string MeshTransientSolver::Solve(
    double t0, double t_end,
    std::function<void(double, const std::vector<std::vector<double>>&)>
        output_cb) {
  // Compute consistent initial conditions; non-fatal if it fails.
  const int ic_flag =
      IDACalcIC(ida_mem_, IDA_YA_YDP_INIT, t0 + opts_.initial_time_step);
  if (ic_flag < 0) {
    std::fprintf(stderr,
                 "[mphys] IDACalcIC failed (flag %d) — proceeding with supplied "
                 "initial conditions\n",
                 ic_flag);
  }

  Scatter(yy_, scratch_y_);
  if (output_cb) output_cb(t0, scratch_y_);

  CheckFlag(IDASetStopTime(ida_mem_, t_end), "IDASetStopTime");

  double t_ret = t0;
  while (t_ret < t_end) {
    const int flag = IDASolve(ida_mem_, t_end, &t_ret, yy_, yp_, IDA_ONE_STEP);
    if (flag < 0) {
      Scatter(yy_, model_.fields());
      return "Solver stopped at t=" + std::to_string(t_ret) +
             " (IDASolve flag " + std::to_string(flag) + ")";
    }
    Scatter(yy_, scratch_y_);
    if (output_cb) output_cb(t_ret, scratch_y_);
    if (flag == IDA_TSTOP_RETURN || t_ret >= t_end) break;
  }

  Scatter(yy_, model_.fields());
  return {};
}

// static
int MeshTransientSolver::ResidualCb(sunrealtype t, N_Vector yy, N_Vector yp,
                                    N_Vector rr, void* user_data) {
  auto& solver = *static_cast<MeshTransientSolver*>(user_data);
  solver.Scatter(yy, solver.scratch_y_);
  solver.Scatter(yp, solver.scratch_ydot_);
  solver.model_.Residual(t, solver.scratch_y_, solver.scratch_ydot_,
                         solver.scratch_rr_);
  solver.Gather(solver.scratch_rr_, rr);
  return 0;
}

void MeshTransientSolver::Scatter(N_Vector nv,
                                  std::vector<std::vector<double>>& fields) const {
  const double* data = N_VGetArrayPointer(nv);
  for (int k = 0; k < n_fields_; ++k) {
    const int base = k * n_cells_;
    for (int i = 0; i < n_cells_; ++i) fields[k][i] = data[base + i];
  }
}

void MeshTransientSolver::Gather(const std::vector<std::vector<double>>& fields,
                                 N_Vector nv) const {
  double* data = N_VGetArrayPointer(nv);
  for (int k = 0; k < n_fields_; ++k) {
    const int base = k * n_cells_;
    for (int i = 0; i < n_cells_; ++i) data[base + i] = fields[k][i];
  }
}

// static
void MeshTransientSolver::CheckFlag(int flag, const char* func_name) {
  if (flag < 0) {
    throw std::runtime_error(std::string(func_name) + " failed with flag " +
                             std::to_string(flag));
  }
}

}  // namespace mphys
