#pragma once

#include <string>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/topology.hpp"

namespace mphys {

// Base class for physics models defined on a face-based Mesh (1D / 2D / 3D).
//
// This is the dimension-independent counterpart to Model: users subclass it,
// declare fields and per-patch boundary conditions in the constructor, and
// implement Residual() using the face-based fvm:: operators. The same model
// runs on any Mesh regardless of dimension.
//
// Solved either to steady state (MeshSteadySolver, Newton-Krylov) or in time
// (MeshTransientSolver, IDA) — the same Residual() serves both. Each field is a
// plain std::vector<double> of length mesh.NCells().
class MeshModel {
 public:
  explicit MeshModel(const Mesh& mesh) : mesh_(mesh) {}
  virtual ~MeshModel() = default;

  // Declare a field (one value per cell). Returns its index, used to index the
  // y / ydot / rr arrays in Residual(). Boundary conditions default to
  // zero-Neumann on every patch until SetBcs() is called.
  int AddField(const std::string& name, double init_value = 0.0);

  // Attach one PatchBc per mesh patch to a field (ordered as Mesh::patches).
  void SetBcs(int field, std::vector<PatchBc> patch_bcs);

  // Mark a field as algebraic (no time derivative) for transient DAE solves.
  // Fields are differential by default.
  void MarkFieldAlgebraic(int field);

  // Form the residual rr = F(t, y, ydot) = 0. y[k] / ydot[k] / rr[k] are the
  // cell values of field k and its time derivative. For a steady solve `ydot`
  // is empty; transient models write e.g. rr[k][c] = ydot[k][c] - rhs[c].
  virtual void Residual(double t,
                        const std::vector<std::vector<double>>& y,
                        const std::vector<std::vector<double>>& ydot,
                        std::vector<std::vector<double>>& rr) = 0;

  const Mesh& mesh() const { return mesh_; }
  int NFields() const { return static_cast<int>(field_names_.size()); }
  const std::string& field_name(int k) const { return field_names_[k]; }
  bool field_is_differential(int k) const { return field_differential_[k]; }
  std::vector<std::vector<double>>& fields() { return fields_; }
  const std::vector<std::vector<double>>& fields() const { return fields_; }
  const std::vector<PatchBc>& bcs(int field) const { return bcs_[field]; }

 protected:
  const Mesh& mesh_;
  std::vector<std::string> field_names_;
  std::vector<std::vector<double>> fields_;  // initial guess, then solution
  std::vector<std::vector<PatchBc>> bcs_;    // [field][patch]
  std::vector<bool> field_differential_;     // [field]; false => algebraic
};

}  // namespace mphys
