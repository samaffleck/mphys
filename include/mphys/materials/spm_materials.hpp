#pragma once

#include <vector>

#include "mphys/materials/database.hpp"
#include "mphys/materials/material.hpp"
#include "mphys/models/spm.hpp"

// ============================================================
// Bridge: materials -> Single Particle Model parameters.
//
// The SPM needs, per electrode, an OCP curve, a maximum concentration, a solid
// diffusivity and a kinetic rate constant — exactly the intrinsic data an
// ElectrodeMaterial carries — plus a nominal electrolyte concentration.  This
// adapter combines two chosen electrode materials and an electrolyte with a
// CellDesign (geometry the materials do not know about) and the operating
// current/temperature to build a fully-populated SpmParameters.
//
// This is the concrete expression of "the model defines the inputs required
// from the material": SpmMaterialSlots() lists what must be selected, and
// MakeSpmParameters() consumes that selection.
// ============================================================

namespace mphys::materials {

// Geometry + initial state that the cell design (not the material) provides.
struct CellDesign {
  double R_n = 5.86e-6;   // negative particle radius          [m]
  double R_p = 5.22e-6;   // positive particle radius          [m]
  double L_n = 85.2e-6;   // negative electrode thickness      [m]
  double L_p = 75.6e-6;   // positive electrode thickness      [m]
  double A = 0.1027;      // electrode cross-sectional area    [m^2]
  double eps_n = 0.75;    // negative active-material fraction [-]
  double eps_p = 0.665;   // positive active-material fraction [-]
  double x_n0 = 0.84;     // initial negative stoichiometry    [-]
  double x_p0 = 0.27;     // initial positive stoichiometry    [-]
};

// Material slots the SPM requires.  A UI can use this to build pickers.
inline std::vector<MaterialSlot> SpmMaterialSlots() {
  return {
      {"Negative electrode", Category::kElectrode, Domain::kNegativeElectrode},
      {"Positive electrode", Category::kElectrode, Domain::kPositiveElectrode},
      {"Electrolyte", Category::kElectrolyte, Domain::kElectrolyte},
  };
}

// Build SpmParameters from a material selection + cell design + operating point.
inline models::SpmParameters MakeSpmParameters(
    const ElectrodeMaterial& negative, const ElectrodeMaterial& positive,
    const ElectrolyteMaterial& electrolyte, const CellDesign& cell,
    double current_A, double temperature_K = 298.15) {
  models::SpmParameters p;

  // Cell geometry.
  p.R_n = cell.R_n;
  p.R_p = cell.R_p;
  p.L_n = cell.L_n;
  p.L_p = cell.L_p;
  p.A = cell.A;
  p.eps_n = cell.eps_n;
  p.eps_p = cell.eps_p;

  // Operating conditions.
  p.I = current_A;
  p.T = temperature_K;

  // Negative electrode material.
  p.cn_max = negative.c_max;
  p.cn0 = cell.x_n0 * negative.c_max;
  p.D_n = negative.diffusivity;
  p.k_n = negative.reaction_rate;
  p.Un = negative.ocp;

  // Positive electrode material.
  p.cp_max = positive.c_max;
  p.cp0 = cell.x_p0 * positive.c_max;
  p.D_p = positive.diffusivity;
  p.k_p = positive.reaction_rate;
  p.Up = positive.ocp;

  // Electrolyte (SPM uses a single nominal concentration).
  p.c_e = electrolyte.c_typical;

  return p;
}

// Convenience overload selecting materials by identifier, e.g.
//   MakeSpmParameters(ElectrodeId::kGraphiteChen2020,
//                     ElectrodeId::kNmc811Chen2020,
//                     ElectrolyteId::kLipf6EcEmcNyman2008, cell, 5.0);
inline models::SpmParameters MakeSpmParameters(
    ElectrodeId negative, ElectrodeId positive, ElectrolyteId electrolyte,
    const CellDesign& cell, double current_A, double temperature_K = 298.15) {
  return MakeSpmParameters(Database::Electrode(negative),
                           Database::Electrode(positive),
                           Database::Electrolyte(electrolyte), cell, current_A,
                           temperature_K);
}

}  // namespace mphys::materials
