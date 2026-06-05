#pragma once

#include <functional>
#include <vector>

#include "mphys/topology.hpp"

namespace mphys {

// Diagonal (Jacobi) preconditioner for the matrix-free Newton-Krylov mesh
// solvers. It extracts the diagonal of the iteration matrix
//   J = dF/dy + cj * dF/dydot
// by finite differences, perturbing whole colour classes of a graph colouring
// of the cell adjacency at once. Because non-adjacent cells don't pollute each
// other's residual, the entire diagonal costs only n_fields * n_colours
// residual evaluations — two colours for a structured 5-/7-point grid.
//
// Preconditioning by 1/diag clusters the eigenvalues of variable-coefficient
// and transient (large cj) systems, cutting GMRES iteration counts. For a
// constant-coefficient uniform operator the diagonal is constant, so it is a
// harmless scalar scaling.
class DiagonalPreconditioner {
 public:
  // Evaluates F (length N, field-major) at a flat state. `ydot` is nullptr for
  // steady solves.
  using EvalFn =
      std::function<void(const double* y, const double* ydot, double* F)>;

  DiagonalPreconditioner(const Mesh& mesh, int n_fields);

  // Recompute the inverse diagonal at base state (y, ydot). `base_F` must be
  // F(y, ydot), which the caller already has. `cj` is IDA's scalar (0 steady).
  void Update(const EvalFn& eval, const double* y, const double* ydot, double cj,
              const double* base_F);

  // z = M^{-1} r (length N; z may alias r).
  void Apply(const double* r, double* z) const;

  int n_colors() const { return n_colors_; }

 private:
  int n_cells_;
  int n_fields_;
  int N_;
  int n_colors_ = 0;
  std::vector<std::vector<int>> cells_by_color_;
  std::vector<double> inv_diag_;  // length N

  // Scratch reused across Update() calls.
  std::vector<double> y_pert_;
  std::vector<double> ydot_pert_;
  std::vector<double> f_pert_;
  std::vector<double> delta_;  // per-cell perturbation, length n_cells_
};

}  // namespace mphys
