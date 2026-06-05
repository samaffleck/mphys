#include "mphys/mesh_model.hpp"

#include <stdexcept>

namespace mphys {

int MeshModel::AddField(const std::string& name, double init_value) {
  const int idx = NFields();
  field_names_.push_back(name);
  fields_.emplace_back(mesh_.NCells(), init_value);
  // Default: zero-Neumann (no flux) on every patch.
  bcs_.emplace_back(mesh_.patches.size(), NeumannBc(0.0));
  field_differential_.push_back(true);
  return idx;
}

void MeshModel::MarkFieldAlgebraic(int field) {
  if (field < 0 || field >= NFields()) {
    throw std::out_of_range("MeshModel::MarkFieldAlgebraic: invalid field index");
  }
  field_differential_[field] = false;
}

void MeshModel::SetBcs(int field, std::vector<PatchBc> patch_bcs) {
  if (field < 0 || field >= NFields()) {
    throw std::out_of_range("MeshModel::SetBcs: invalid field index");
  }
  if (patch_bcs.size() != mesh_.patches.size()) {
    throw std::invalid_argument(
        "MeshModel::SetBcs: expected one PatchBc per mesh patch");
  }
  bcs_[field] = std::move(patch_bcs);
}

}  // namespace mphys
