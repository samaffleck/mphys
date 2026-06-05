#pragma once

#include <string>
#include <vector>

namespace mphys {

// Declarative description of the parameters and boundary conditions a model
// exposes. The GUI renders its configuration panel from this schema instead of
// hardcoding per-model widgets, and ModelBuildContext uses it to expand
// per-domain scalars to the per-cell / per-face data each model needs.

// Where a parameter's value is defined.
enum class ParamScope {
  kPerDomain,  // one value per geometry domain (e.g. diffusivity)
  kGlobal,     // a single value shared across the whole model (reserved)
};

// How a per-domain parameter is interpolated to mesh faces when the model
// needs a face-centred field (length n_cells + 1).
enum class FaceInterp {
  kNone,        // cell-centred only
  kHarmonic,    // resistance-weighted harmonic mean (diffusivities)
  kArithmetic,  // distance-weighted arithmetic mean (velocities)
};

struct ParamSpec {
  std::string key;            // stable identifier, e.g. "D"
  std::string label;          // UI label, e.g. "Diffusivity"
  std::string unit;           // UI unit (UTF-8 ok), e.g. "m\xc2\xb2/s"
  double      default_value = 0.0;
  ParamScope  scope = ParamScope::kPerDomain;
  FaceInterp  face_interp = FaceInterp::kNone;
};

// One selectable option for a boundary slot — e.g. "Dirichlet"/"Neumann", or
// "Pressure"/"Velocity". The model's factory interprets the chosen option.
struct BcOption {
  std::string label;          // e.g. "Dirichlet"
  std::string unit;           // unit of the value field for this option
};

// A boundary the model exposes (typically the "left"/"right" terminal nodes).
struct BcSlot {
  std::string key;            // "left" or "right"
  std::string label;          // "Boundary 1"
  std::vector<BcOption> options;
  int    default_option = 0;
  double default_value = 0.0;
};

struct ParamSchema {
  std::vector<ParamSpec> params;
  std::vector<BcSlot>    boundaries;
};

}  // namespace mphys
