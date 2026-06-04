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
// Steady-state only for now (solved by MeshSteadySolver via Newton-Krylov).
// Each field is a plain std::vector<double> of length mesh.NCells().
class MeshModel {
 public:
  explicit MeshModel(const Mesh& mesh) : mesh_(mesh) {}
  virtual ~MeshModel() = default;

  // Declare a field (one value per cell). Returns its index, used to index the
  // y / rr arrays in Residual(). Boundary conditions default to zero-Neumann on
  // every patch until SetBcs() is called.
  int AddField(const std::string& name, double init_value = 0.0);

  // Attach one PatchBc per mesh patch to a field (ordered as Mesh::patches).
  void SetBcs(int field, std::vector<PatchBc> patch_bcs);

  // Form the steady residual rr = F(y) = 0. y[k] / rr[k] are the cell values of
  // field k. Use mesh() and bcs(k) with the fvm:: face operators.
  virtual void Residual(const std::vector<std::vector<double>>& y,
                        std::vector<std::vector<double>>& rr) = 0;

  const Mesh& mesh() const { return mesh_; }
  int NFields() const { return static_cast<int>(field_names_.size()); }
  const std::string& field_name(int k) const { return field_names_[k]; }
  std::vector<std::vector<double>>& fields() { return fields_; }
  const std::vector<std::vector<double>>& fields() const { return fields_; }
  const std::vector<PatchBc>& bcs(int field) const { return bcs_[field]; }

 protected:
  const Mesh& mesh_;
  std::vector<std::string> field_names_;
  std::vector<std::vector<double>> fields_;  // initial guess, then solution
  std::vector<std::vector<PatchBc>> bcs_;    // [field][patch]
};

}  // namespace mphys
