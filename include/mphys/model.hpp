#pragma once

#include <string>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/mesh.hpp"
#include "mphys/state_vector.hpp"

namespace mphys {

// Base class for user-defined physics models.
//
// Users subclass Model, declare fields and boundary conditions in the
// constructor, and implement Residual() to express the governing equations
// using fvm:: operators.
class Model {
 public:
  explicit Model(const Mesh1D& mesh, StateVector& sv);
  virtual ~Model() = default;

  // Declare a spatial field (one value per cell).
  // Returns the field index for use in Residual() to index into y, ydot, rr.
  int AddField(const std::string& name, double init_value = 0.0);

  // Declare a 0D algebraic variable.
  // Returns the algebraic index for use in Residual() to index into y_alg, rr_alg.
  int AddAlgebraic(const std::string& name, double init_value = 0.0);

  // Associate boundary conditions with a field, identified by index or name.
  void SetBcs(int field_index, FieldBcs bcs);
  void SetBcs(const std::string& field_name, FieldBcs bcs);

  // Evaluate the residual F(t, y, ydot) = 0.
  // For transient problems: all arrays are populated.
  // For steady-state (KINSOL): ydot and ydot_alg are empty vectors.
  virtual void Residual(double t,
                        const std::vector<Field>& y,
                        const std::vector<Field>& ydot,
                        const std::vector<double>& y_alg,
                        const std::vector<double>& ydot_alg,
                        std::vector<Field>& rr,
                        std::vector<double>& rr_alg) = 0;

  // Accessors used by solvers.
  const Mesh1D& mesh() const { return mesh_; }
  StateVector& state_vector() { return sv_; }
  std::vector<Field>& fields() { return fields_; }
  std::vector<double>& algebraics() { return algebraics_; }
  const std::vector<FieldBcs>& bcs() const { return bcs_; }

 protected:
  const Mesh1D& mesh_;
  StateVector& sv_;
  std::vector<Field> fields_;      // indexed same as sv_.var_names_
  std::vector<double> algebraics_; // indexed same as sv_.algebraic_names_
  std::vector<FieldBcs> bcs_;      // indexed same as fields_
};

}  // namespace mphys
