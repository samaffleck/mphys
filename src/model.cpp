#include "mphys/model.hpp"

#include <stdexcept>

namespace mphys {

Model::Model(const Mesh1D& mesh, StateVector& sv) : mesh_(mesh), sv_(sv) {}

int Model::AddField(const std::string& name, double init_value) {
  sv_.AddField(name);
  fields_.emplace_back(name, mesh_.n_cells, init_value);
  bcs_.emplace_back();  // default: Neumann(0) on both sides
  return static_cast<int>(fields_.size()) - 1;
}

int Model::AddAlgebraic(const std::string& name, double init_value) {
  sv_.AddAlgebraic(name);
  algebraics_.push_back(init_value);
  return static_cast<int>(algebraics_.size()) - 1;
}

void Model::SetBcs(int field_index, FieldBcs bcs) {
  if (field_index < 0 || field_index >= static_cast<int>(bcs_.size())) {
    throw std::out_of_range("SetBcs: field index out of range");
  }
  bcs_[field_index] = bcs;
}

void Model::SetBcs(const std::string& field_name, FieldBcs bcs) {
  for (int k = 0; k < static_cast<int>(fields_.size()); ++k) {
    if (fields_[k].name == field_name) {
      bcs_[k] = bcs;
      return;
    }
  }
  throw std::invalid_argument("SetBcs: field '" + field_name + "' not registered");
}

}  // namespace mphys
