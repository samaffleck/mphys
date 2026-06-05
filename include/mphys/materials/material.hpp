#pragma once

#include <functional>
#include <string>

// ============================================================
// Materials database — core types.
//
// A *material* describes the intrinsic, geometry-independent physical
// properties of a battery component (e.g. the equilibrium potential of NMC811,
// or the conductivity of 1 M LiPF6 in EC:EMC).  Materials are grouped by
// Category (electrode / electrolyte) and carry a "typical" Domain hint so the
// GUI can offer only sensible choices for a given slot.
//
// Materials deliberately do *not* contain cell-design parameters (electrode
// thickness, particle radius, area, active-material fraction) or operating
// conditions (current, temperature).  Those belong to the cell/experiment.  A
// physics model declares which material slots it needs (see MaterialSlot and
// each model's RequiredMaterials()) and an adapter combines the chosen
// materials with a cell design to produce the model's parameter struct.
// ============================================================

namespace mphys::materials {

// Broad grouping used to populate material pickers.
enum class Category {
  kElectrode,
  kElectrolyte,
};

// Where a material physically lives in a cell.  Electrode materials carry a
// "typical" domain (graphite -> negative, NMC/LFP -> positive) used only as a
// soft filter for the UI; the same active material may in principle be used in
// either electrode.
enum class Domain {
  kNegativeElectrode,
  kPositiveElectrode,
  kSeparator,
  kElectrolyte,
  kAny,
};

// Property that depends on a single scalar — e.g. an electrode open-circuit
// potential U(x), where x = c_surf / c_max is the surface stoichiometry [-].
using ScalarProperty = std::function<double(double)>;

// Property that depends on local state (salt concentration c [mol/m^3]) and
// temperature (T [K]) — e.g. electrolyte conductivity kappa(c, T) [S/m].
using StateProperty = std::function<double(double c, double T)>;

// Wrap a constant value as a property (handy for properties whose dependence we
// don't yet model).
inline ScalarProperty Constant(double value) {
  return [value](double) { return value; };
}
inline StateProperty Constant2(double value) {
  return [value](double, double) { return value; };
}

// ------------------------------------------------------------
// Electrode active material (graphite, NMC811, LFP, ...).
//
// Properties are expressed against the local stoichiometry x = c / c_max in
// [0, 1].  Defaults of zero mark "not provided"; a model that needs a property
// should validate it.
// ------------------------------------------------------------
struct ElectrodeMaterial {
  std::string name;                     // display name, e.g. "NMC811 (Chen 2020)"
  std::string chemistry;                // short tag, e.g. "NMC811", "graphite"
  Domain typical_domain = Domain::kAny;
  std::string reference;                // literature source for the parameters

  // --- Thermodynamics ---
  ScalarProperty ocp;                   // open-circuit potential U(x) [V vs Li/Li+]
  ScalarProperty entropic_coeff = Constant(0.0);  // dU/dT(x) [V/K]
  double c_max = 0.0;                   // max Li concentration [mol/m^3]
  double x_min = 0.0;                   // usable stoichiometry window [-]
  double x_max = 1.0;

  // --- Solid-phase transport ---
  double diffusivity = 0.0;             // D_s [m^2/s] (reference temperature)
  double activation_energy_D = 0.0;     // Arrhenius E_a for D_s [J/mol] (0 = none)

  // --- Kinetics (Butler-Volmer) ---
  double reaction_rate = 0.0;           // rate constant k (see model for units)
  double activation_energy_k = 0.0;     // Arrhenius E_a for k [J/mol] (0 = none)

  // --- Bulk physical properties ---
  double density = 0.0;                 // [kg/m^3]
  double molar_mass = 0.0;              // [kg/mol]
};

// ------------------------------------------------------------
// Electrolyte (salt + solvent), e.g. 1 M LiPF6 in EC:EMC.
//
// Transport properties are inherently functions of salt concentration and
// temperature, so they are StateProperty callables.
// ------------------------------------------------------------
struct ElectrolyteMaterial {
  std::string name;                     // display name
  std::string reference;                // literature source

  StateProperty conductivity;           // kappa(c, T) [S/m]
  StateProperty diffusivity;            // D_e(c, T) [m^2/s]
  StateProperty transference;           // cation transference number t_+(c, T) [-]
  StateProperty thermodynamic_factor = Constant2(1.0);  // 1 + dln f_+/dln c [-]
  double c_typical = 1000.0;            // nominal salt concentration [mol/m^3]
};

// ------------------------------------------------------------
// A material requirement declared by a physics model.  A model exposes a list
// of these (see e.g. SpmMaterialSlots()) so a UI can build the right pickers
// and validate a selection before a run.
// ------------------------------------------------------------
struct MaterialSlot {
  std::string role;       // human label, e.g. "Negative electrode"
  Category category;      // which kind of material fits this slot
  Domain domain;          // compatible domain (used to filter choices)
};

}  // namespace mphys::materials
