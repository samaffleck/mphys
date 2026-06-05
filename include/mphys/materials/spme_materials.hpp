#pragma once

#include <vector>

#include "mphys/materials/database.hpp"
#include "mphys/materials/material.hpp"
#include "mphys/materials/spm_materials.hpp"
#include "mphys/models/spme.hpp"

// ============================================================
// Bridge: materials -> Single Particle Model with Electrolyte parameters.
//
// SPMe reuses the SPM particle/voltage "core" (built by MakeSpmParameters) and
// adds a macroscopic electrolyte transport equation.  The leading-order SPMe is
// a *constant-property* electrolyte model: it needs scalar values for the salt
// diffusivity D_e, ionic conductivity kappa_e and transference number t_+.  An
// ElectrolyteMaterial carries those as functions of concentration/temperature,
// so this adapter samples them once at the nominal salt concentration and the
// operating temperature.  The separator thickness, electrolyte porosities,
// Bruggeman exponent and solid-phase conductivities are cell design, supplied
// via SpmeCellDesign.
// ============================================================

namespace mphys::materials {

// Cell design for SPMe: the SPM particle geometry plus the electrolyte-side
// geometry that the materials do not know about.
struct SpmeCellDesign {
  CellDesign particles;     // particle radii, thicknesses, area, eps, x0 (SPM)

  double L_s = 12e-6;       // separator thickness                      [m]
  double eps_e_n = 0.25;    // negative electrode electrolyte fraction  [-]
  double eps_e_s = 0.47;    // separator electrolyte fraction           [-]
  double eps_e_p = 0.335;   // positive electrode electrolyte fraction  [-]
  double brugg = 1.5;       // Bruggeman exponent                       [-]
  double sigma_n = 215.0;   // negative solid-phase conductivity        [S/m]
  double sigma_p = 0.18;    // positive solid-phase conductivity        [S/m]
};

// Material slots the SPMe requires (same shape as SPM, but the electrolyte's
// transport functions are actually consumed here).
inline std::vector<MaterialSlot> SpmeMaterialSlots() {
  return {
      {"Negative electrode", Category::kElectrode, Domain::kNegativeElectrode},
      {"Positive electrode", Category::kElectrode, Domain::kPositiveElectrode},
      {"Electrolyte", Category::kElectrolyte, Domain::kElectrolyte},
  };
}

// Build SpmeParameters from a material selection + cell design + operating point.
inline models::SpmeParameters MakeSpmeParameters(
    const ElectrodeMaterial& negative, const ElectrodeMaterial& positive,
    const ElectrolyteMaterial& electrolyte, const SpmeCellDesign& cell,
    double current_A, double temperature_K = 298.15) {
  models::SpmeParameters p;

  // Particle/voltage core reuses the SPM bridge (OCPs, c_max, D_s, kinetics,
  // geometry, operating point, nominal electrolyte concentration).
  p.core = MakeSpmParameters(negative, positive, electrolyte, cell.particles,
                             current_A, temperature_K);

  // Electrolyte-side geometry (cell design).
  p.L_s = cell.L_s;
  p.eps_e_n = cell.eps_e_n;
  p.eps_e_s = cell.eps_e_s;
  p.eps_e_p = cell.eps_e_p;
  p.brugg = cell.brugg;
  p.sigma_n = cell.sigma_n;
  p.sigma_p = cell.sigma_p;

  // Electrolyte transport: sample the material's functions at the nominal salt
  // concentration and operating temperature (leading-order constant-property).
  const double ce = electrolyte.c_typical;
  p.ce0 = ce;
  p.D_e = electrolyte.diffusivity(ce, temperature_K);
  p.kappa_e = electrolyte.conductivity(ce, temperature_K);
  p.t_plus = electrolyte.transference(ce, temperature_K);

  return p;
}

// Convenience overload selecting materials by identifier.
inline models::SpmeParameters MakeSpmeParameters(
    ElectrodeId negative, ElectrodeId positive, ElectrolyteId electrolyte,
    const SpmeCellDesign& cell, double current_A,
    double temperature_K = 298.15) {
  return MakeSpmeParameters(Database::Electrode(negative),
                            Database::Electrode(positive),
                            Database::Electrolyte(electrolyte), cell, current_A,
                            temperature_K);
}

}  // namespace mphys::materials
