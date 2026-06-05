#pragma once

#include <functional>
#include <optional>
#include <string>
#include <vector>

#include <nvector/nvector_serial.h>

#include "mphys/diagonal_preconditioner.hpp"
#include "mphys/mesh_model.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/sundials_types.hpp"

namespace mphys {

// Transient solver for a MeshModel using IDA with a matrix-free Krylov linear
// solver (SPGMR). Like MeshSteadySolver it forms no Jacobian matrix, so it runs
// on 1D/2D/3D meshes without a banded-bandwidth assumption.
//
// The model's Residual(t, y, ydot, rr) is integrated as a DAE F(t,y,ydot)=0.
// Fields marked algebraic (MeshModel::MarkFieldAlgebraic) get id = 0 so IDA can
// compute consistent initial conditions. Unknowns are laid out field-major:
// index = field * n_cells + cell.
class MeshTransientSolver {
 public:
  MeshTransientSolver(MeshModel& model, const SolverOptions& opts,
                      SUNContext sunctx);
  ~MeshTransientSolver() = default;

  MeshTransientSolver(const MeshTransientSolver&) = delete;
  MeshTransientSolver& operator=(const MeshTransientSolver&) = delete;

  // Integrate from t0 to t_end. The callback (if given) fires with the field
  // values after the initial state and each accepted step. Returns an empty
  // string on success, or a message if the integrator stopped early.
  std::string Solve(
      double t0, double t_end,
      std::function<void(double t, const std::vector<std::vector<double>>&)>
          output_cb = {});

  // Total GMRES iterations across the integration (diagnostic / benchmarking).
  long NumLinearIterations() const;

 private:
  static int ResidualCb(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector rr,
                        void* user_data);
  static int PrecSetupCb(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr,
                         sunrealtype cj, void* user_data);
  static int PrecSolveCb(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr,
                         N_Vector rvec, N_Vector zvec, sunrealtype cj,
                         sunrealtype delta, void* user_data);
  static void CheckFlag(int flag, const char* func_name);

  void Scatter(N_Vector nv, std::vector<std::vector<double>>& fields) const;
  void Gather(const std::vector<std::vector<double>>& fields, N_Vector nv) const;

  MeshModel& model_;
  SolverOptions opts_;
  int n_cells_;
  int n_fields_;
  int N_;

  IdaMem ida_mem_;
  SunVector yy_;
  SunVector yp_;
  SunVector id_;
  SunLinearSolver ls_;

  std::optional<DiagonalPreconditioner> precond_;

  std::vector<std::vector<double>> scratch_y_;
  std::vector<std::vector<double>> scratch_ydot_;
  std::vector<std::vector<double>> scratch_rr_;
};

}  // namespace mphys
