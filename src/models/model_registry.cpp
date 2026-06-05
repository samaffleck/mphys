#include "mphys/models/model_registry.hpp"

#include <algorithm>
#include <utility>

namespace mphys {

void ModelRegistry::Register(ModelInfo info) {
  models_.push_back(std::move(info));
}

std::vector<std::string> ModelRegistry::Packages() const {
  std::vector<std::string> packages;
  for (const auto& m : models_)
    if (std::find(packages.begin(), packages.end(), m.package) == packages.end())
      packages.push_back(m.package);
  return packages;
}

std::vector<const ModelInfo*> ModelRegistry::InPackage(
    const std::string& package) const {
  std::vector<const ModelInfo*> out;
  for (const auto& m : models_)
    if (m.package == package) out.push_back(&m);
  return out;
}

const ModelInfo* ModelRegistry::Find(const std::string& id) const {
  for (const auto& m : models_)
    if (m.id == id) return &m;
  return nullptr;
}

ModelRegistry& BuiltinModels() {
  static ModelRegistry registry = [] {
    ModelRegistry r;
    RegisterBuiltinModels(r);
    return r;
  }();
  return registry;
}

}  // namespace mphys
