#include "mphys/state_vector.hpp"

#include <cassert>
#include <stdexcept>

#include <nvector/nvector_serial.h>

namespace mphys {

StateVector::StateVector(int n_cells) : n_cells_(n_cells) {}

void StateVector::AddField(std::string name) {
  var_names_.push_back(std::move(name));
  field_is_differential_.push_back(true);
  ++n_vars_;
}

void StateVector::AddAlgebraic(std::string name) {
  algebraic_names_.push_back(std::move(name));
  ++n_algebraics_;
}

void StateVector::MarkFieldAlgebraic(const std::string& name) {
  for (int k = 0; k < n_vars_; ++k) {
    if (var_names_[k] == name) {
      field_is_differential_[k] = false;
      return;
    }
  }
  throw std::invalid_argument("MarkFieldAlgebraic: field '" + name + "' not registered");
}

void StateVector::Scatter(N_Vector nv, std::vector<Field>& fields,
                          std::vector<double>& algebraics) const {
  assert(static_cast<int>(fields.size()) == n_vars_);
  assert(static_cast<int>(algebraics.size()) == n_algebraics_);

  const double* data = N_VGetArrayPointer(nv);
  for (int i = 0; i < n_cells_; ++i) {
    for (int k = 0; k < n_vars_; ++k) {
      fields[k][i] = data[i * n_vars_ + k];
    }
  }
  const int spatial_size = n_cells_ * n_vars_;
  for (int j = 0; j < n_algebraics_; ++j) {
    algebraics[j] = data[spatial_size + j];
  }
}

void StateVector::Gather(const std::vector<Field>& fields,
                         const std::vector<double>& algebraics, N_Vector nv) const {
  assert(static_cast<int>(fields.size()) == n_vars_);
  assert(static_cast<int>(algebraics.size()) == n_algebraics_);

  double* data = N_VGetArrayPointer(nv);
  for (int i = 0; i < n_cells_; ++i) {
    for (int k = 0; k < n_vars_; ++k) {
      data[i * n_vars_ + k] = fields[k][i];
    }
  }
  const int spatial_size = n_cells_ * n_vars_;
  for (int j = 0; j < n_algebraics_; ++j) {
    data[spatial_size + j] = algebraics[j];
  }
}

void StateVector::FillIdVector(N_Vector id) const {
  double* data = N_VGetArrayPointer(id);
  // Spatial block
  for (int i = 0; i < n_cells_; ++i) {
    for (int k = 0; k < n_vars_; ++k) {
      data[i * n_vars_ + k] = field_is_differential_[k] ? 1.0 : 0.0;
    }
  }
  // Algebraic block: all 0.0
  const int spatial_size = n_cells_ * n_vars_;
  for (int j = 0; j < n_algebraics_; ++j) {
    data[spatial_size + j] = 0.0;
  }
}

void StateVector::AllocateScratch(std::vector<Field>& fields,
                                  std::vector<double>& algebraics) const {
  fields.resize(n_vars_);
  for (int k = 0; k < n_vars_; ++k) {
    fields[k] = Field(var_names_[k], n_cells_, 0.0);
  }
  algebraics.assign(n_algebraics_, 0.0);
}

}  // namespace mphys
