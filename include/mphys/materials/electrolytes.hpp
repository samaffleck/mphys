#pragma once

#include <cmath>

#include "mphys/materials/material.hpp"

// ============================================================
// Electrolyte data.
//
// Transport properties are functions of salt concentration c [mol/m^3] and
// temperature T [K].  The polynomial fits below are quoted for concentration in
// mol/L, so each builder rescales c internally (ce = c / 1000).  Temperature
// dependence is not yet modelled for these sets; the T argument is accepted for
// forward compatibility (Arrhenius corrections can be added later).
// ============================================================

namespace mphys::materials {

// 1 M LiPF6 in EC:EMC (3:7 wt), Nyman et al. 2008 — the standard PyBaMM
// electrolyte for the LG M50 / Chen 2020 parameter set.
inline ElectrolyteMaterial Lipf6_EcEmc_Nyman2008() {
  ElectrolyteMaterial m;
  m.name = "LiPF6 in EC:EMC (Nyman 2008)";
  m.reference = "Nyman et al., Electrochim. Acta 53 6356 (2008)";
  m.conductivity = [](double c, double) {
    const double ce = c / 1000.0;  // mol/L
    return 0.1297 * std::pow(ce, 3) - 2.51 * std::pow(ce, 1.5) + 3.329 * ce;
  };
  m.diffusivity = [](double c, double) {
    const double ce = c / 1000.0;  // mol/L
    return (8.794e-11 * ce * ce - 3.972e-10 * ce + 4.862e-10);
  };
  m.transference = Constant2(0.2594);
  m.thermodynamic_factor = Constant2(1.0);
  m.c_typical = 1000.0;
  return m;
}

// 1 M LiPF6 in EC:DMC, Marquis et al. 2019 (Capiglia-based diffusivity).
inline ElectrolyteMaterial Lipf6_EcDmc_Marquis2019() {
  ElectrolyteMaterial m;
  m.name = "LiPF6 in EC:DMC (Marquis 2019)";
  m.reference = "Marquis et al., J. Electrochem. Soc. 166 A3693 (2019)";
  m.conductivity = [](double c, double) {
    const double ce = c / 1000.0;  // mol/L
    return 0.0911 + 1.9101 * ce - 1.052 * ce * ce + 0.1554 * std::pow(ce, 3);
  };
  m.diffusivity = [](double c, double) {
    const double ce = c / 1000.0;  // mol/L
    return 5.34e-10 * std::exp(-0.65 * ce);
  };
  m.transference = Constant2(0.4);
  m.thermodynamic_factor = Constant2(1.0);
  m.c_typical = 1000.0;
  return m;
}

}  // namespace mphys::materials
