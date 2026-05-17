#pragma once

#include <functional>
#include <vector>

#include <sundials/sundials_types.h>

#include "mphys/field.hpp"
#include "mphys/model.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sundials_types.hpp"

namespace mphys {

// Wraps SUNDIALS IDA to integrate a Model forward in time.
//
// Solve() drives the time loop and fires an optional callback after each
// accepted step so the caller can record or post-process results.
class TransientSolver {
 public:
  TransientSolver(Model& model, const SolverOptions& opts, SUNContext sunctx);
  // Default destructor — all SUNDIALS memory freed by RAII members.
  ~TransientSolver() = default;

  TransientSolver(const TransientSolver&) = delete;
  TransientSolver& operator=(const TransientSolver&) = delete;

  // Returns an empty string on success, or a warning/error message if the
  // solver stopped early (e.g. convergence failure).  Partial results are
  // still delivered via output_cb for every step that succeeded.
  std::string Solve(double t0, double t_end,
                    std::function<void(double t, const std::vector<Field>&,
                                       const std::vector<double>&)>
                        output_cb = {});

 private:
  static int ResidualCb(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector rr,
                        void* user_data);

  static void CheckFlag(int flag, const char* func_name);

  Model& model_;
  SolverOptions opts_;

  IdaMem          ida_mem_;
  SunMatrix       A_;
  SunLinearSolver ls_;
  SunVector       yy_;
  SunVector       yp_;
  SunVector       id_;

  // Scratch buffers allocated once to avoid heap allocation in the hot path.
  std::vector<Field>  scratch_y_;
  std::vector<Field>  scratch_ydot_;
  std::vector<Field>  scratch_rr_;
  std::vector<double> scratch_alg_;
  std::vector<double> scratch_aldot_;
  std::vector<double> scratch_ralg_;
};

}  // namespace mphys
