#pragma once

#include <string>
#include <vector>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>

#include "mphys/field.hpp"

namespace mphys {

// Maps named Field objects to positions in a flat SUNDIALS N_Vector.
//
// Layout: [var0_cell0, var1_cell0, ..., varK_cell0,
//          var0_cell1, var1_cell1, ..., varK_cell1, ...
//          alg0, alg1, ...]
//
// Half-bandwidth of the banded Jacobian = n_vars (number of spatial fields).
class StateVector {
 public:
  explicit StateVector(int n_cells);

  // Register a spatial field (one value per cell). Must be called before Scatter/Gather.
  void AddField(std::string name);

  // Register a 0D algebraic variable appended after all spatial data.
  void AddAlgebraic(std::string name);

  // Mark all cells of a spatial field as algebraic (id = 0) for IDASetId.
  void MarkFieldAlgebraic(const std::string& name);

  int TotalSize() const { return n_cells_ * n_vars_ + n_algebraics_; }
  int NCells() const { return n_cells_; }
  int NVars() const { return n_vars_; }
  int NAlgebraics() const { return n_algebraics_; }

  // Copy data from an N_Vector into Field/algebraic arrays.
  void Scatter(N_Vector nv, std::vector<Field>& fields,
               std::vector<double>& algebraics) const;

  // Copy data from Field/algebraic arrays into an N_Vector.
  void Gather(const std::vector<Field>& fields,
              const std::vector<double>& algebraics, N_Vector nv) const;

  // Fill an existing N_Vector with 1.0 (differential) or 0.0 (algebraic) per component.
  void FillIdVector(N_Vector id) const;

  // Allocate scratch Field/algebraic arrays sized to match registered variables.
  void AllocateScratch(std::vector<Field>& fields,
                       std::vector<double>& algebraics) const;

 private:
  int n_cells_;
  int n_vars_ = 0;
  int n_algebraics_ = 0;
  std::vector<std::string> var_names_;
  std::vector<std::string> algebraic_names_;
  // id_flags_[k] = true if spatial field k is differential (default), false if algebraic.
  std::vector<bool> field_is_differential_;
};

}  // namespace mphys
