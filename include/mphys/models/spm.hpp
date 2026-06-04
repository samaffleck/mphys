#pragma once

#include <cmath>
#include <functional>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/fvm_operators.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/state_vector.hpp"

// ============================================================
// Single Particle Model (SPM)
//
// Reference: PyBaMM SPM notebook
//   https://docs.pybamm.org/en/stable/source/examples/notebooks/models/SPM.html
//
// Each electrode is represented by a single spherical particle in which
// lithium diffuses (Fick's law in spherical coordinates):
//
//   dc_s,k/dt = (1/r²) d/dr( r² D_s,k dc_s,k/dr )    k in {n, p}
//
// Boundary conditions:
//   dc_s,k/dr = 0                       at r = 0   (symmetry / no flux)
//   -D_s,k dc_s,k/dr = ∓ I/(F a_k L_k)  at r = R_k (surface, sign per electrode)
//
// where a_k = 3 ε_k / R_k is the specific surface area, L_k the electrode
// thickness, ε_k the active-material volume fraction, A the electrode area,
// and I the applied current (discharge positive).
//
// To keep both particles on a single shared mesh, each particle is solved on a
// *normalised* radial coordinate r = r̃/R_k ∈ [0,1].  The diffusion coefficient
// then scales to D_k/R_k², and the surface gradient (the Neumann BC) becomes
//
//   dc_n/dr|_{r=1} = -R_n² i_app / (3 D_n F ε_n L_n)   (negative electrode)
//   dc_p/dr|_{r=1} = +R_p² i_app / (3 D_p F ε_p L_p)   (positive electrode)
//
// with i_app = I/A.  These signs make the negative particle deplete and the
// positive particle fill during discharge, and reproduce exactly the
// couliometric balance  d⟨c_k⟩/dt = ∓ i_app / (F ε_k L_k).
//
// The terminal voltage is an algebraic output:
//
//   V = U_p(y_surf) - U_n(x_surf)
//       - (2RT/F) asinh( i_app / (2 j0_p a_p L_p) )
//       - (2RT/F) asinh( i_app / (2 j0_n a_n L_n) )
//
// with stoichiometries x = c/c_max and exchange-current densities
//   j0_k = F k_k sqrt(c_e) sqrt(c_surf) sqrt(c_max - c_surf).
// ============================================================

namespace mphys::models {

// Default open-circuit potentials (Chen et al. 2020, LG M50 graphite / NMC811).
// Argument is the surface stoichiometry x = c_surf / c_max.
inline double DefaultOcpNegative(double x) {
  return 1.9793 * std::exp(-39.3631 * x) + 0.2482
       - 0.0909 * std::tanh(29.8538 * (x - 0.1234))
       - 0.04478 * std::tanh(14.9159 * (x - 0.2769))
       - 0.0205 * std::tanh(30.4444 * (x - 0.6103));
}

inline double DefaultOcpPositive(double y) {
  return -0.8090 * y + 4.4875
       - 0.0428 * std::tanh(18.5138 * (y - 0.5542))
       - 17.7326 * std::tanh(15.7890 * (y - 0.3117))
       + 17.5842 * std::tanh(15.9308 * (y - 0.3120));
}

// All SPM parameters.  Defaults are representative LG M50 (Chen 2020) values.
struct SpmParameters {
  // Physical constants
  static constexpr double F  = 96485.33212;   // Faraday constant   [C/mol]
  static constexpr double Rg = 8.314462618;   // gas constant       [J/mol/K]

  // Geometry
  double R_n = 5.86e-6;   // negative particle radius            [m]
  double R_p = 5.22e-6;   // positive particle radius            [m]
  double L_n = 85.2e-6;   // negative electrode thickness        [m]
  double L_p = 75.6e-6;   // positive electrode thickness        [m]
  double A   = 0.1027;    // electrode cross-sectional area      [m²]
  double eps_n = 0.75;    // negative active-material fraction   [-]
  double eps_p = 0.665;   // positive active-material fraction   [-]

  // Solid-phase transport
  double D_n = 3.3e-14;   // negative solid diffusivity          [m²/s]
  double D_p = 4.0e-15;   // positive solid diffusivity          [m²/s]

  // Concentrations
  double cn_max = 33133.0;            // negative max concentration   [mol/m³]
  double cp_max = 63104.0;            // positive max concentration   [mol/m³]
  double cn0    = 0.84 * 33133.0;     // initial negative concentration
  double cp0    = 0.27 * 63104.0;     // initial positive concentration
  double c_e    = 1000.0;             // electrolyte concentration    [mol/m³]

  // Kinetics (Butler-Volmer rate constants)
  double k_n = 6.48e-7;   // negative reaction rate constant
  double k_p = 3.42e-6;   // positive reaction rate constant

  // Operating conditions
  double I = 5.0;         // applied current, discharge positive [A]
  double T = 298.15;      // temperature                         [K]

  // Open-circuit potentials (overridable).
  std::function<double(double)> Un = DefaultOcpNegative;
  std::function<double(double)> Up = DefaultOcpPositive;

  // Derived quantities
  double a_n()   const { return 3.0 * eps_n / R_n; }   // specific area [1/m]
  double a_p()   const { return 3.0 * eps_p / R_p; }
  double i_app() const { return I / A; }               // current density [A/m²]
  double Dn_eff() const { return D_n / (R_n * R_n); }  // on r ∈ [0,1]
  double Dp_eff() const { return D_p / (R_p * R_p); }
};

// Exchange-current density [A/m²] for one electrode.
inline double SpmExchangeCurrent(double k, double c_e, double c_surf,
                                 double c_max) {
  const double a = std::max(c_surf, 0.0);
  const double b = std::max(c_max - c_surf, 0.0);
  return SpmParameters::F * k * std::sqrt(c_e) * std::sqrt(a) * std::sqrt(b);
}

// Terminal voltage [V] from the two surface concentrations.
inline double SpmVoltage(const SpmParameters& p, double cn_surf, double cp_surf) {
  const double x = cn_surf / p.cn_max;
  const double y = cp_surf / p.cp_max;
  const double j0n = SpmExchangeCurrent(p.k_n, p.c_e, cn_surf, p.cn_max);
  const double j0p = SpmExchangeCurrent(p.k_p, p.c_e, cp_surf, p.cp_max);
  const double rt2f = 2.0 * SpmParameters::Rg * p.T / SpmParameters::F;
  const double i = p.i_app();

  const double eta_n = rt2f * std::asinh(i / (2.0 * j0n * p.a_n() * p.L_n));
  const double eta_p = rt2f * std::asinh(i / (2.0 * j0p * p.a_p() * p.L_p));

  return p.Up(y) - p.Un(x) - eta_p - eta_n;
}

// SPM as an mphys Model: two diffusion fields ("c_n", "c_p") on a normalised
// spherical mesh r ∈ [0,1], plus an algebraic terminal-voltage output ("V").
//
// Construct the mesh with MakeUniformMesh1D(0.0, 1.0, n, CoordSystem::kSpherical).
class SpmModel : public Model {
 public:
  SpmModel(const Mesh1D& mesh, StateVector& sv, SpmParameters params)
      : Model(mesh, sv), p_(params) {
    cn_ = AddField("c_n", p_.cn0);
    cp_ = AddField("c_p", p_.cp0);
    V_  = AddAlgebraic("V", SpmVoltage(p_, p_.cn0, p_.cp0));

    // Surface gradients (Neumann BC at r = 1) — constant for constant current.
    gn_ = -p_.R_n * p_.R_n * p_.i_app() /
          (3.0 * p_.D_n * SpmParameters::F * p_.eps_n * p_.L_n);
    gp_ = +p_.R_p * p_.R_p * p_.i_app() /
          (3.0 * p_.D_p * SpmParameters::F * p_.eps_p * p_.L_p);

    SetBcs(cn_, {NeumannBc(0.0), NeumannBc(gn_)});
    SetBcs(cp_, {NeumannBc(0.0), NeumannBc(gp_)});
  }

  void Residual(double,
                const std::vector<Field>& y,
                const std::vector<Field>& ydot,
                const std::vector<double>& y_alg,
                const std::vector<double>&,
                std::vector<Field>& rr,
                std::vector<double>& rr_alg) override {
    rr[cn_] = fvm::Ddt(ydot[cn_])
            - p_.Dn_eff() * fvm::Laplacian(y[cn_], 1.0, mesh_, bcs_[cn_]);
    rr[cp_] = fvm::Ddt(ydot[cp_])
            - p_.Dp_eff() * fvm::Laplacian(y[cp_], 1.0, mesh_, bcs_[cp_]);

    const double cn_s = SurfaceValue(y[cn_], gn_);
    const double cp_s = SurfaceValue(y[cp_], gp_);
    rr_alg[V_] = y_alg[V_] - SpmVoltage(p_, cn_s, cp_s);
  }

  // Extrapolate a field to the particle surface (r = 1) using its Neumann
  // gradient: c_surf = c[last] + g * dx/2.
  double SurfaceValue(const Field& c, double g) const {
    const int n = c.NCells();
    return c[n - 1] + g * 0.5 * mesh_.dx[n - 1];
  }

  const SpmParameters& params() const { return p_; }
  int c_n_index() const { return cn_; }
  int c_p_index() const { return cp_; }
  int voltage_index() const { return V_; }

 private:
  SpmParameters p_;
  int cn_ = 0, cp_ = 1, V_ = 0;
  double gn_ = 0.0, gp_ = 0.0;
};

}  // namespace mphys::models
