#pragma once

#include <string>
#include <vector>

#include "mphys/models/model_info.hpp"

namespace mphys {

// Catalogue of available models, grouped by package. Replaces the hardcoded
// physics enum + dropdown in the GUI: the GUI builds its model tree by querying
// the registry, and constructs a model via the matching ModelInfo::factory.
class ModelRegistry {
 public:
  void Register(ModelInfo info);

  // All registered models, in registration order.
  const std::vector<ModelInfo>& All() const { return models_; }

  // Distinct package names, in first-seen order.
  std::vector<std::string> Packages() const;

  // Models belonging to a package, in registration order.
  std::vector<const ModelInfo*> InPackage(const std::string& package) const;

  // Lookup by stable id; nullptr if not found.
  const ModelInfo* Find(const std::string& id) const;

 private:
  std::vector<ModelInfo> models_;
};

// Registers all built-in models (Transport, Fluid Flow, ...).
void RegisterBuiltinModels(ModelRegistry& registry);

// Process-wide registry, populated with the built-in models on first access.
ModelRegistry& BuiltinModels();

}  // namespace mphys
