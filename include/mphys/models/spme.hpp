#pragma once

#include <cmath>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/fvm_operators.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/models/spm.hpp"
#include "mphys/state_vector.hpp"

// ============================================================
// Single Particle Model with Electrolyte (SPMe)
//
// Reference: PyBaMM SPMe notebook
//   https://docs.pybamm.org/en/stable/source/examples/notebooks/models/SPMe.html
//   (Marquis et al. 2019, "An asymptotic derivation of a single particle model
//    with electrolyte", J. Electrochem. Soc. 166 A3693.)
//
// The SPMe keeps the two SPM spherical-particle diffusion equations unchanged
// and adds a macroscopic electrolyte-concentration equation across the cell
// sandwich (negative electrode | separator | positive electrode) together with
// the leading-order electrolyte corrections to the terminal voltage.  Because
// the asymptotic reduction decouples the particle and electrolyte transport,
// the two physics never couple inside the residual — they only meet in the
// scalar voltage output.  This lets us solve them on two independent meshes
// that merely share the same cell count:
//
//   * particles  c_n, c_p  — normalised spherical mesh r ∈ [0,1]   (as in SPM)
//   * electrolyte c_e      — Cartesian mesh x ∈ [0, L_n+L_s+L_p]
//
// Electrolyte transport (porous-electrode mass balance):
//
//   ε(x) dc_e/dt = d/dx( D_e^eff(x) dc_e/dx ) + S(x)
//
//   D_e^eff(x) = ε(x)^b D_e          (Bruggeman correction)
//   S(x) = +(1-t+) i_app / (F L_n)   in the negative electrode
//        = 0                          in the separator
//        = -(1-t+) i_app / (F L_p)   in the positive electrode
//
// with no-flux boundaries  dc_e/dx = 0  at x = 0 and x = L.  The net source
// integrates to zero, so the porosity-weighted average ⟨ε c_e⟩ is conserved.
//
// Terminal voltage (leading-order SPMe):
//
//   V = U_p(y) - U_n(x)                              (open-circuit)
//     - (2RT/F) asinh( i_app / (2 j0_p a_p L_p) )    (positive reaction)
//     - (2RT/F) asinh( i_app / (2 j0_n a_n L_n) )    (negative reaction)
//     + 2(1-t+)(RT/F) ( ln c̄_e,p - ln c̄_e,n )       (concentration overpot.)
//     - i_app ( L_n/(3κ_n) + L_s/κ_s + L_p/(3κ_p) )  (electrolyte ohmic)
//     - (i_app/3)( L_n/σ_n + L_p/σ_p )               (solid ohmic)
//
// where the exchange-current densities j0_k now use the electrode-averaged
// electrolyte concentration c̄_e,k, κ_k = ε_e,k^b κ_e is the effective
// electrolyte conductivity, and σ_k the solid-phase conductivity.  Setting the
// electrolyte corrections to zero recovers the SPM voltage exactly.
// ============================================================

namespace mphys::models {

// Electrolyte-specific parameters layered on top of the SPM particle/voltage
// core.  Defaults follow the LG M50 (Chen 2020) parameter set used by SPM.
struct SpmeParameters {
  SpmParameters core;     // particle geometry, transport, kinetics, operating

  // Separator
  double L_s = 12e-6;     // separator thickness                    [m]

  // Electrolyte volume fractions (porosities) per region
  double eps_e_n = 0.25;  // negative electrode electrolyte fraction [-]
  double eps_e_s = 0.47;  // separator electrolyte fraction          [-]
  double eps_e_p = 0.335; // positive electrode electrolyte fraction [-]

  // Electrolyte transport / thermodynamics
  double D_e     = 1.769e-10;  // electrolyte salt diffusivity       [m²/s]
  double kappa_e = 1.0;        // electrolyte ionic conductivity     [S/m]
  double t_plus  = 0.2594;     // cation transference number         [-]
  double brugg   = 1.5;        // Bruggeman exponent                 [-]

  // Solid-phase conductivities
  double sigma_n = 215.0;      // negative solid conductivity        [S/m]
  double sigma_p = 0.18;       // positive solid conductivity        [S/m]

  // Initial / reference electrolyte concentration
  double ce0 = 1000.0;         // initial electrolyte concentration  [mol/m³]

  // Derived per-region effective transport (Bruggeman).
  double De_n() const { return std::pow(eps_e_n, brugg) * D_e; }
  double De_s() const { return std::pow(eps_e_s, brugg) * D_e; }
  double De_p() const { return std::pow(eps_e_p, brugg) * D_e; }
  double kappa_n() const { return std::pow(eps_e_n, brugg) * kappa_e; }
  double kappa_s() const { return std::pow(eps_e_s, brugg) * kappa_e; }
  double kappa_p() const { return std::pow(eps_e_p, brugg) * kappa_e; }

  // Total sandwich length.
  double L() const { return core.L_n + L_s + core.L_p; }
};

// A three-region Cartesian electrolyte mesh whose region interfaces fall
// exactly on cell faces, plus a per-cell region label (0=neg, 1=sep, 2=pos).
struct ElectrolyteMesh {
  Mesh1D mesh;
  std::vector<int> region;   // length n_cells
};

// Build the electrolyte mesh from uniform sub-meshes for each region so that
// porosity and source terms are exactly piecewise-constant.
inline ElectrolyteMesh MakeSpmeElectrolyteMesh(const SpmeParameters& p,
                                               int n_n, int n_s, int n_p) {
  const double Ln = p.core.L_n, Ls = p.L_s, Lp = p.core.L_p;
  const int n = n_n + n_s + n_p;

  ElectrolyteMesh em;
  Mesh1D& m = em.mesh;
  m.n_cells = n;
  m.coord_system = CoordSystem::kCartesian;
  m.face_positions.resize(n + 1);
  m.cell_centres.resize(n);
  m.dx.resize(n);
  em.region.resize(n);

  const double dn = Ln / n_n, ds = Ls / n_s, dp = Lp / n_p;
  m.face_positions[0] = 0.0;
  for (int i = 0; i < n; ++i) {
    double w;
    if (i < n_n)            { w = dn; em.region[i] = 0; }
    else if (i < n_n + n_s) { w = ds; em.region[i] = 1; }
    else                    { w = dp; em.region[i] = 2; }
    m.face_positions[i + 1] = m.face_positions[i] + w;
    m.cell_centres[i] = 0.5 * (m.face_positions[i] + m.face_positions[i + 1]);
    m.dx[i] = w;
  }
  m.df.resize(n - 1);
  for (int i = 0; i < n - 1; ++i) {
    m.df[i] = m.cell_centres[i + 1] - m.cell_centres[i];
  }
  return em;
}

// Terminal voltage [V] including the electrolyte corrections.  c̄_e,n / c̄_e,p
// are the electrode-averaged electrolyte concentrations.
inline double SpmeVoltage(const SpmeParameters& p, double cn_surf,
                          double cp_surf, double ce_n_avg, double ce_p_avg) {
  const SpmParameters& c = p.core;
  const double F = SpmParameters::F;
  const double rt_f = SpmParameters::Rg * c.T / F;
  const double i = c.i_app();

  // Open-circuit potentials.
  const double x = cn_surf / c.cn_max;
  const double y = cp_surf / c.cp_max;
  const double ocv = c.Up(y) - c.Un(x);

  // Reaction overpotentials — exchange current uses the local electrolyte
  // concentration in each electrode.
  const double j0n = SpmExchangeCurrent(c.k_n, ce_n_avg, cn_surf, c.cn_max);
  const double j0p = SpmExchangeCurrent(c.k_p, ce_p_avg, cp_surf, c.cp_max);
  const double eta_n = 2.0 * rt_f * std::asinh(i / (2.0 * j0n * c.a_n() * c.L_n));
  const double eta_p = 2.0 * rt_f * std::asinh(i / (2.0 * j0p * c.a_p() * c.L_p));

  // Concentration overpotential (assumes unit thermodynamic factor).
  const double eta_c = 2.0 * (1.0 - p.t_plus) * rt_f *
                       (std::log(ce_p_avg) - std::log(ce_n_avg));

  // Electrolyte ohmic drop (leading-order SPMe).
  const double dphi_e = -i * (c.L_n / (3.0 * p.kappa_n()) + p.L_s / p.kappa_s() +
                              c.L_p / (3.0 * p.kappa_p()));

  // Solid-phase ohmic drop.
  const double dphi_s = -(i / 3.0) * (c.L_n / p.sigma_n + c.L_p / p.sigma_p);

  return ocv - eta_p - eta_n + eta_c + dphi_e + dphi_s;
}

// SPMe as an mphys Model.  Fields: "c_n", "c_p" (spherical particle mesh) and
// "c_e" (Cartesian electrolyte mesh); algebraic output "V".  All three fields
// share the same cell count.
//
// Construct with:
//   auto em   = MakeSpmeElectrolyteMesh(p, n_n, n_s, n_p);
//   auto mesh = MakeUniformMesh1D(0,1, em.mesh.n_cells, CoordSystem::kSpherical);
//   StateVector sv(em.mesh.n_cells);
//   SpmeModel model(mesh, em, sv, p);
class SpmeModel : public Model {
 public:
  SpmeModel(const Mesh1D& particle_mesh, ElectrolyteMesh electrolyte_mesh,
            StateVector& sv, SpmeParameters params)
      : Model(particle_mesh, sv), p_(params), em_(std::move(electrolyte_mesh)) {
    const SpmParameters& c = p_.core;

    cn_ = AddField("c_n", c.cn0);
    cp_ = AddField("c_p", c.cp0);
    ce_ = AddField("c_e", p_.ce0);
    V_  = AddAlgebraic("V", 0.0);

    // Particle surface gradients (Neumann BC at r = 1) — as in SPM.
    gn_ = -c.R_n * c.R_n * c.i_app() /
          (3.0 * c.D_n * SpmParameters::F * c.eps_n * c.L_n);
    gp_ = +c.R_p * c.R_p * c.i_app() /
          (3.0 * c.D_p * SpmParameters::F * c.eps_p * c.L_p);
    SetBcs(cn_, {NeumannBc(0.0), NeumannBc(gn_)});
    SetBcs(cp_, {NeumannBc(0.0), NeumannBc(gp_)});

    // Electrolyte: no-flux at both current collectors.
    SetBcs(ce_, {NeumannBc(0.0), NeumannBc(0.0)});

    BuildElectrolyteFields();

    // Seed the algebraic voltage with a consistent initial value.
    algebraics_[V_] = SpmeVoltage(p_, c.cn0, c.cp0, p_.ce0, p_.ce0);
  }

  void Residual(double,
                const std::vector<Field>& y,
                const std::vector<Field>& ydot,
                const std::vector<double>& y_alg,
                const std::vector<double>&,
                std::vector<Field>& rr,
                std::vector<double>& rr_alg) override {
    const SpmParameters& c = p_.core;

    // Particle diffusion on the shared normalised spherical mesh.
    rr[cn_] = fvm::Ddt(ydot[cn_]) -
              c.Dn_eff() * fvm::Laplacian(y[cn_], 1.0, mesh_, bcs_[cn_]);
    rr[cp_] = fvm::Ddt(ydot[cp_]) -
              c.Dp_eff() * fvm::Laplacian(y[cp_], 1.0, mesh_, bcs_[cp_]);

    // Electrolyte diffusion on the Cartesian sandwich mesh:
    //   ε dc_e/dt = div(D_e^eff grad c_e) + S
    Field diff = fvm::Laplacian(y[ce_], De_face_, em_.mesh, bcs_[ce_]);
    for (int i = 0; i < em_.mesh.n_cells; ++i) {  // rr[ce_] is pre-sized scratch
      rr[ce_][i] = eps_e_[i] * ydot[ce_][i] - diff[i] - source_[i];
    }

    // Voltage output.
    const double cn_s = SurfaceValue(y[cn_], gn_);
    const double cp_s = SurfaceValue(y[cp_], gp_);
    const double ce_n = RegionAverage(y[ce_], 0);
    const double ce_p = RegionAverage(y[ce_], 2);
    rr_alg[V_] = y_alg[V_] - SpmeVoltage(p_, cn_s, cp_s, ce_n, ce_p);
  }

  // Extrapolate a particle field to its surface (r = 1) via its Neumann grad.
  double SurfaceValue(const Field& c, double g) const {
    const int n = c.NCells();
    return c[n - 1] + g * 0.5 * mesh_.dx[n - 1];
  }

  // Volume (dx-weighted) average of an electrolyte field over one region.
  double RegionAverage(const Field& ce, int region) const {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < em_.mesh.n_cells; ++i) {
      if (em_.region[i] != region) continue;
      num += em_.mesh.dx[i] * ce[i];
      den += em_.mesh.dx[i];
    }
    return den > 0.0 ? num / den : p_.ce0;
  }

  const SpmeParameters& params() const { return p_; }
  const ElectrolyteMesh& electrolyte_mesh() const { return em_; }
  int c_n_index() const { return cn_; }
  int c_p_index() const { return cp_; }
  int c_e_index() const { return ce_; }
  int voltage_index() const { return V_; }

 private:
  // Precompute per-cell porosity & source and per-face effective diffusivity.
  void BuildElectrolyteFields() {
    const SpmParameters& c = p_.core;
    const int n = em_.mesh.n_cells;

    eps_e_.assign(n, 0.0);
    source_.assign(n, 0.0);
    std::vector<double> De_cell(n, 0.0);

    const double s_scale = (1.0 - p_.t_plus) * c.i_app() / SpmParameters::F;
    for (int i = 0; i < n; ++i) {
      switch (em_.region[i]) {
        case 0:  // negative electrode: Li+ produced on discharge
          eps_e_[i] = p_.eps_e_n; De_cell[i] = p_.De_n();
          source_[i] = +s_scale / c.L_n; break;
        case 1:  // separator
          eps_e_[i] = p_.eps_e_s; De_cell[i] = p_.De_s();
          source_[i] = 0.0; break;
        default: // positive electrode: Li+ consumed
          eps_e_[i] = p_.eps_e_p; De_cell[i] = p_.De_p();
          source_[i] = -s_scale / c.L_p; break;
      }
    }

    // Face diffusivity via harmonic mean of neighbouring cells (correct flux
    // continuity across the discontinuous region interfaces).  Boundary faces
    // are unused (zero-flux BC sets the boundary gradient directly).
    De_face_ = Field("De_face", n + 1, 0.0);
    De_face_[0] = De_cell[0];
    for (int i = 1; i < n; ++i) {
      De_face_[i] = 2.0 * De_cell[i - 1] * De_cell[i] /
                    (De_cell[i - 1] + De_cell[i]);
    }
    De_face_[n] = De_cell[n - 1];
  }

  SpmeParameters p_;
  ElectrolyteMesh em_;
  int cn_ = 0, cp_ = 1, ce_ = 2, V_ = 0;
  double gn_ = 0.0, gp_ = 0.0;

  std::vector<double> eps_e_;   // per-cell porosity
  std::vector<double> source_;  // per-cell source term
  Field De_face_;               // per-face effective diffusivity
};

}  // namespace mphys::models
