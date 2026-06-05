#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/models/param_schema.hpp"
#include "mphys/state_vector.hpp"

namespace mphys {

// Which solver a model is integrated with.
enum class SolverKind { kTransient, kSteady };

// A user's chosen configuration for one boundary slot.
struct BcChoice {
  int    option = 0;   // index into the slot's BcSlot::options
  double value  = 0.0;
};

// Everything a model factory needs to build a Model from a configured geometry
// and parameter set. Per-domain scalar parameters are expanded to per-cell or
// per-face data through the helper methods, centralising the interpolation that
// each model would otherwise repeat.
class ModelBuildContext {
 public:
  using ParamLookup = std::function<double(int domain, const std::string& key)>;
  using BcLookup    = std::function<BcChoice(const std::string& slot)>;

  ModelBuildContext(const Mesh1D& mesh, StateVector& sv,
                    std::vector<int> cell_domain,
                    ParamLookup param_lookup, BcLookup bc_lookup);

  const Mesh1D& mesh() const { return mesh_; }
  StateVector&  sv()   const { return sv_; }
  int n_domains() const;

  // Raw per-domain value of a parameter.
  double DomainParam(int domain, const std::string& key) const;

  // Per-domain scalar -> per-cell vector (piecewise constant).
  std::vector<double> CellParam(const std::string& key) const;

  // Per-domain scalar -> face field (length n_cells + 1) using the given interp.
  Field FaceParam(const std::string& key, FaceInterp interp) const;

  // The user's choice for a boundary slot ("left"/"right").
  BcChoice Bc(const std::string& slot) const;

 private:
  const Mesh1D& mesh_;
  StateVector&  sv_;
  std::vector<int> cell_domain_;
  ParamLookup param_lookup_;
  BcLookup    bc_lookup_;
};

// Metadata + factory for one runnable model. Registered with a ModelRegistry.
struct ModelInfo {
  std::string id;            // stable key, e.g. "transport.conv_diff_reaction"
  std::string package;       // e.g. "Transport"
  std::string group;         // optional sub-label within a package, "" if none
  std::string name;          // e.g. "Convection-Diffusion-Reaction"
  std::string description;
  SolverKind  solver = SolverKind::kTransient;
  ParamSchema schema;
  std::function<std::unique_ptr<Model>(ModelBuildContext&)> factory;

  // Models whose geometry, parameters or results don't fit the generic
  // schema-driven path set this; the GUI then dispatches to bespoke panels and
  // a custom runner keyed by id (the factory above is left unset). The Single
  // Particle Model is the first such model.
  bool custom_gui = false;
};

}  // namespace mphys
