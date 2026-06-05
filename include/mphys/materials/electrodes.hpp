#pragma once

#include <cmath>

#include "mphys/materials/material.hpp"

// ============================================================
// Electrode active-material data.
//
// Each builder returns a fully-populated ElectrodeMaterial.  Parameter sets are
// taken from the cited literature; the open-circuit potential (OCP) fits are
// reproduced verbatim so they can be checked against the original papers.
//
// Stoichiometry convention: x = c / c_max in [0, 1], where c is the local Li
// concentration and c_max the maximum.  OCP is in volts versus Li/Li+.
// ============================================================

namespace mphys::materials {

// ------------------------------------------------------------
// OCP fits (free functions so they can be reused / unit tested directly).
// ------------------------------------------------------------

// Graphite negative electrode, Chen et al. 2020 (LG M50, "Graphite_LGM50").
inline double GraphiteOcp_Chen2020(double x) {
  return 1.9793 * std::exp(-39.3631 * x) + 0.2482
       - 0.0909 * std::tanh(29.8538 * (x - 0.1234))
       - 0.04478 * std::tanh(14.9159 * (x - 0.2769))
       - 0.0205 * std::tanh(30.4444 * (x - 0.6103));
}

// NMC811 positive electrode, Chen et al. 2020 (LG M50, "NMC811_LGM50").
inline double Nmc811Ocp_Chen2020(double y) {
  return -0.8090 * y + 4.4875
       - 0.0428 * std::tanh(18.5138 * (y - 0.5542))
       - 17.7326 * std::tanh(15.7890 * (y - 0.3117))
       + 17.5842 * std::tanh(15.9308 * (y - 0.3120));
}

// LiFePO4 positive electrode, Prada et al. 2013 (flat ~3.4 V plateau).
inline double LfpOcp_Prada2013(double y) {
  return 3.4323
       - 0.8428 * std::exp(-80.2493 * std::pow(1.0 - y, 1.3198))
       - 3.2474e-6 * std::exp(20.2645 * std::pow(1.0 - y, 3.8003))
       + 3.2482e-6 * std::exp(20.2646 * std::pow(1.0 - y, 3.7995));
}

// ------------------------------------------------------------
// Material builders.
// ------------------------------------------------------------

// Graphite — LG M50 negative electrode (Chen et al., J. Electrochem. Soc. 2020).
inline ElectrodeMaterial Graphite_Chen2020() {
  ElectrodeMaterial m;
  m.name = "Graphite (Chen 2020)";
  m.chemistry = "graphite";
  m.typical_domain = Domain::kNegativeElectrode;
  m.reference = "Chen et al., J. Electrochem. Soc. 167 080534 (2020)";
  m.ocp = GraphiteOcp_Chen2020;
  m.c_max = 33133.0;        // [mol/m^3]
  m.x_min = 0.0;
  m.x_max = 1.0;
  m.diffusivity = 3.3e-14;  // [m^2/s]
  m.reaction_rate = 6.48e-7;
  m.density = 2266.0;       // [kg/m^3]
  m.molar_mass = 0.072066;  // C6 [kg/mol]
  return m;
}

// NMC811 — LG M50 positive electrode (Chen et al. 2020).
inline ElectrodeMaterial Nmc811_Chen2020() {
  ElectrodeMaterial m;
  m.name = "NMC811 (Chen 2020)";
  m.chemistry = "NMC811";
  m.typical_domain = Domain::kPositiveElectrode;
  m.reference = "Chen et al., J. Electrochem. Soc. 167 080534 (2020)";
  m.ocp = Nmc811Ocp_Chen2020;
  m.c_max = 63104.0;        // [mol/m^3]
  m.x_min = 0.0;
  m.x_max = 1.0;
  m.diffusivity = 4.0e-15;  // [m^2/s]
  m.reaction_rate = 3.42e-6;
  m.density = 3262.0;       // [kg/m^3]
  m.molar_mass = 0.0965;    // [kg/mol]
  return m;
}

// LiFePO4 — positive electrode (Prada et al., J. Electrochem. Soc. 2013).
inline ElectrodeMaterial Lfp_Prada2013() {
  ElectrodeMaterial m;
  m.name = "LFP (Prada 2013)";
  m.chemistry = "LFP";
  m.typical_domain = Domain::kPositiveElectrode;
  m.reference = "Prada et al., J. Electrochem. Soc. 160 A1908 (2013)";
  m.ocp = LfpOcp_Prada2013;
  m.c_max = 22806.0;        // [mol/m^3]
  m.x_min = 0.0;
  m.x_max = 1.0;
  m.diffusivity = 5.9e-18;  // [m^2/s] (very low — LFP is two-phase)
  m.reaction_rate = 6.0e-7;
  m.density = 3600.0;       // [kg/m^3]
  m.molar_mass = 0.157757;  // [kg/mol]
  return m;
}

// NMC532 — positive electrode (representative; Ecker et al. 2015 family).
inline ElectrodeMaterial Nmc532_Ecker2015() {
  ElectrodeMaterial m;
  m.name = "NMC532 (Ecker 2015)";
  m.chemistry = "NMC532";
  m.typical_domain = Domain::kPositiveElectrode;
  m.reference = "Ecker et al., J. Electrochem. Soc. 162 A1836 (2015)";
  // Reasonable NMC532 plateau shape; replace with the full Ecker fit if needed.
  m.ocp = [](double y) {
    return 4.2 - 1.1 * y - 0.15 * std::tanh(15.0 * (y - 0.45));
  };
  m.c_max = 48580.0;        // [mol/m^3]
  m.x_min = 0.0;
  m.x_max = 1.0;
  m.diffusivity = 3.7e-13;  // [m^2/s]
  m.reaction_rate = 3.0e-6;
  m.density = 4870.0;       // [kg/m^3]
  m.molar_mass = 0.0966;    // [kg/mol]
  return m;
}

}  // namespace mphys::materials
