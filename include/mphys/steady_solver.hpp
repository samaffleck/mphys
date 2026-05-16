#pragma once

#include <vector>

#include <sundials/sundials_types.h>

#include "mphys/field.hpp"
#include "mphys/model.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sundials_types.hpp"

namespace mphys {

// Newton globalisation strategy for SteadySolver.
//   kPure       — pure Newton (robust for linear and mildly nonlinear problems)
//   kLineSearch — Newton with line search (better for strongly nonlinear problems)
enum class NewtonStrategy { kPure, kLineSearch };

// Wraps SUNDIALS KINSOL to find the steady-state solution of a Model.
//
// Calls model.Residual() with empty ydot arrays so the same physics
// implementation works for both transient and steady-state solves.
// The solution is written back into model.fields() on success.
class SteadySolver {
 public:
  SteadySolver(Model& model, const SolverOptions& opts, SUNContext sunctx);
  ~SteadySolver() = default;

  SteadySolver(const SteadySolver&) = delete;
  SteadySolver& operator=(const SteadySolver&) = delete;

  void Solve(NewtonStrategy strategy = NewtonStrategy::kPure);

 private:
  static int SystemCb(N_Vector uu, N_Vector fval, void* user_data);

  static void CheckFlag(int flag, const char* func_name);

  Model& model_;
  SolverOptions opts_;

  KinMem          kin_mem_;
  SunMatrix       A_;
  SunLinearSolver ls_;
  SunVector       u_;

  // Scratch buffers for SystemCb (allocated once).
  std::vector<Field>  scratch_y_;
  std::vector<Field>  scratch_rr_;
  std::vector<double> scratch_alg_;
  std::vector<double> scratch_ralg_;
  // Empty ydot arrays passed to Residual() to signal steady-state.
  std::vector<Field>  empty_ydot_;
  std::vector<double> empty_aldot_;
};

}  // namespace mphys
