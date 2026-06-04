#pragma once

#include <utility>
#include <vector>

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

// Boundary condition applied over a whole boundary patch (which may contain many
// faces in 2D/3D). The condition is either uniform across the patch or carries
// one value per face. Used by the face-based fvm:: operators, indexed by
// Face::patch with the per-face value selected by Face::patch_face.
//
// Implicitly constructible from a BoundaryCondition so that a uniform patch can
// be written simply as DirichletBc(v) / NeumannBc(v).
struct PatchBc {
  BcType type = BcType::kNeumann;
  double uniform = 0.0;
  std::vector<double> values;  // empty => `uniform` applies to every face

  PatchBc() = default;
  PatchBc(BoundaryCondition bc) : type(bc.type), uniform(bc.value) {}  // implicit
  PatchBc(BcType t, std::vector<double> v) : type(t), values(std::move(v)) {}

  double Value(int patch_face) const {
    return values.empty() ? uniform : values[patch_face];
  }
};

}  // namespace mphys
