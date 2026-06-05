#include "mphys/diagonal_preconditioner.hpp"

#include <cmath>

namespace mphys {

DiagonalPreconditioner::DiagonalPreconditioner(const Mesh& mesh, int n_fields)
    : n_cells_(mesh.NCells()),
      n_fields_(n_fields),
      N_(n_cells_ * n_fields),
      inv_diag_(N_, 1.0),
      y_pert_(N_, 0.0),
      ydot_pert_(N_, 0.0),
      f_pert_(N_, 0.0),
      delta_(n_cells_, 0.0) {
  // Build cell adjacency from internal faces.
  std::vector<std::vector<int>> adj(n_cells_);
  for (const Face& f : mesh.faces) {
    if (f.neighbour >= 0) {
      adj[f.owner].push_back(f.neighbour);
      adj[f.neighbour].push_back(f.owner);
    }
  }

  // Greedy graph colouring: adjacent cells get different colours, so each
  // colour class is an independent set in the residual stencil.
  std::vector<int> color(n_cells_, -1);
  std::vector<char> used;
  for (int c = 0; c < n_cells_; ++c) {
    used.assign(n_colors_ + 1, 0);
    for (int nb : adj[c]) {
      if (color[nb] >= 0 && color[nb] < static_cast<int>(used.size())) {
        used[color[nb]] = 1;
      }
    }
    int col = 0;
    while (col < static_cast<int>(used.size()) && used[col]) ++col;
    color[c] = col;
    if (col + 1 > n_colors_) n_colors_ = col + 1;
  }

  cells_by_color_.assign(n_colors_, {});
  for (int c = 0; c < n_cells_; ++c) cells_by_color_[color[c]].push_back(c);
}

void DiagonalPreconditioner::Update(const EvalFn& eval, const double* y,
                                    const double* ydot, double cj,
                                    const double* base_F) {
  constexpr double kEps = 1e-7;

  std::copy(y, y + N_, y_pert_.begin());
  const bool transient = ydot != nullptr;
  if (transient) std::copy(ydot, ydot + N_, ydot_pert_.begin());

  for (int k = 0; k < n_fields_; ++k) {
    const int base = k * n_cells_;
    for (const auto& group : cells_by_color_) {
      // Perturb every cell of this colour in field k simultaneously.
      for (int c : group) {
        const int idx = base + c;
        const double d = kEps * (std::abs(y[idx]) + kEps);
        delta_[c] = d;
        y_pert_[idx] += d;
        if (transient) ydot_pert_[idx] += cj * d;
      }

      eval(y_pert_.data(), transient ? ydot_pert_.data() : nullptr,
           f_pert_.data());

      // Read the diagonal entry for each perturbed cell and restore the state.
      for (int c : group) {
        const int idx = base + c;
        const double diag = (f_pert_[idx] - base_F[idx]) / delta_[c];
        inv_diag_[idx] = (diag != 0.0) ? 1.0 / diag : 1.0;
        y_pert_[idx] = y[idx];
        if (transient) ydot_pert_[idx] = ydot[idx];
      }
    }
  }
}

void DiagonalPreconditioner::Apply(const double* r, double* z) const {
  for (int i = 0; i < N_; ++i) z[i] = inv_diag_[i] * r[i];
}

}  // namespace mphys
