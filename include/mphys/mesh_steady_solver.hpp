#pragma once

#include <vector>

#include <nvector/nvector_serial.h>

#include "mphys/mesh_model.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sundials_types.hpp"

namespace mphys {

// Steady-state solver for a MeshModel using a matrix-free Newton-Krylov method
// (KINSOL + SPGMR). No Jacobian matrix is formed: Jacobian-vector products come
// from finite-difference directional derivatives, so the cost and memory scale
// with the mesh rather than with a dense/banded bandwidth. This is what lets
// the same solver handle 1D/2D/3D meshes; a sparse-direct solver (KLU) can be
// dropped in later without touching the model layer.
//
// The unknowns are laid out field-major: index = field * n_cells + cell.
class MeshSteadySolver {
 public:
  MeshSteadySolver(MeshModel& model, const SolverOptions& opts, SUNContext sunctx);
  ~MeshSteadySolver() = default;

  MeshSteadySolver(const MeshSteadySolver&) = delete;
  MeshSteadySolver& operator=(const MeshSteadySolver&) = delete;

  // Run the nonlinear solve. On success the solution is written back into
  // model.fields(). Throws on failure.
  void Solve();

 private:
  static int SystemCb(N_Vector uu, N_Vector fval, void* user_data);
  static void CheckFlag(int flag, const char* func_name);

  void Scatter(N_Vector nv, std::vector<std::vector<double>>& fields) const;
  void Gather(const std::vector<std::vector<double>>& fields, N_Vector nv) const;

  MeshModel& model_;
  SolverOptions opts_;
  int n_cells_;
  int n_fields_;
  int N_;

  KinMem kin_mem_;
  SunVector u_;
  SunLinearSolver ls_;

  std::vector<std::vector<double>> scratch_y_;
  std::vector<std::vector<double>> scratch_rr_;
};

}  // namespace mphys
