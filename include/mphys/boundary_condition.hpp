#pragma once

namespace mphys {

enum class BcType { kDirichlet, kNeumann };

struct BoundaryCondition {
  BcType type = BcType::kNeumann;
  double value = 0.0;
};

struct FieldBcs {
  BoundaryCondition left;
  BoundaryCondition right;
};

// Convenience constructors
inline BoundaryCondition DirichletBc(double value) {
  return {BcType::kDirichlet, value};
}

inline BoundaryCondition NeumannBc(double value) {
  return {BcType::kNeumann, value};
}

}  // namespace mphys
