/**
 * mphys_gui.cpp  —  COMSOL-inspired desktop front-end for mphys
 *
 * Layout:
 *   left  (~18%)  Model Builder tree
 *   centre (~30%) Configuration panel
 *   right  (~52%) Geometry View (always visible, CAD canvas)
 */

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include <GLFW/glfw3.h>
#include "themes.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/html5.h>
#else
#include "tinyfiledialogs.h"
#endif

#include "mphys/boundary_condition.hpp"
#include "mphys/fvm_operators.hpp"
#include "mphys/materials/database.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/models/spm.hpp"
#include "mphys/models/spme.hpp"
#include "mphys/sim_result.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/steady_solver.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/transient_solver.hpp"

// ============================================================
// Expression evaluator — supports +,-,*,/,^,(,), unary minus,
// constants pi/e, and common math functions.
// ============================================================

namespace expr {
  struct Parser {
    const char* p;
    void skip() { while (*p == ' ' || *p == '\t') ++p; }

    double primary() {
      skip();
      if (*p == '-') { ++p; return -primary(); }
      if (*p == '+') { ++p; return  primary(); }
      if (*p == '(') {
        ++p;
        double v = expr_val();
        skip(); if (*p == ')') ++p;
        return v;
      }
      if (std::isalpha((unsigned char)*p)) {
        const char* s = p;
        while (std::isalnum((unsigned char)*p) || *p == '_') ++p;
        std::string name(s, p);
        skip();
        if (*p == '(') {
          ++p;
          double a = expr_val();
          skip(); if (*p == ')') ++p;
          if (name=="sqrt")            return std::sqrt(a);
          if (name=="sin")             return std::sin(a);
          if (name=="cos")             return std::cos(a);
          if (name=="tan")             return std::tan(a);
          if (name=="log"||name=="ln") return std::log(a);
          if (name=="log10")           return std::log10(a);
          if (name=="exp")             return std::exp(a);
          if (name=="abs")             return std::abs(a);
          if (name=="floor")           return std::floor(a);
          if (name=="ceil")            return std::ceil(a);
          return std::numeric_limits<double>::quiet_NaN();
        }
        if (name=="pi"||name=="PI") return 3.14159265358979323846;
        if (name=="e" ||name=="E")  return 2.71828182845904523536;
        return std::numeric_limits<double>::quiet_NaN();
      }
      char* end;
      double v = std::strtod(p, &end);
      if (end == p) return std::numeric_limits<double>::quiet_NaN();
      p = end;
      return v;
    }

    double power() {
      double b = primary();
      skip();
      if (*p == '^') { ++p; return std::pow(b, power()); }
      return b;
    }

    double term() {
      double v = power();
      for (;;) {
        skip();
        if      (*p == '*') { ++p; v *= power(); }
        else if (*p == '/') { ++p; double d = power(); v = d != 0.0 ? v/d : std::numeric_limits<double>::quiet_NaN(); }
        else break;
      }
      return v;
    }

    double expr_val() {
      double v = term();
      for (;;) {
        skip();
        if      (*p == '+') { ++p; v += term(); }
        else if (*p == '-') { ++p; v -= term(); }
        else break;
      }
      return v;
    }
  };

  inline double eval(const char* s) {
    if (!s || !*s) return 0.0;
    Parser pr{s};
    return pr.expr_val();
  }
} // namespace expr

// ============================================================
// Per-widget expression input state
// ============================================================

struct ExprState { char buf[128] = ""; bool editing = false; };
static std::unordered_map<ImGuiID, ExprState> g_expr_states;

// Expression-aware float input. When idle shows formatted number; while focused
// lets the user type a full expression. On commit, evaluates and updates *v.
// If the field is cleared to empty, *cleared is set to true (value unchanged).
static bool ExprInputFloat(const char* id_str, float* v, bool* cleared = nullptr,
                            const char* fmt = "%.6g") {
  ImGuiID wid = ImGui::GetID(id_str);
  ExprState& st = g_expr_states[wid];

  if (!st.editing)
    snprintf(st.buf, sizeof(st.buf), fmt, (double)*v);

  ImGui::InputText(id_str, st.buf, sizeof(st.buf));

  if (ImGui::IsItemActivated())           st.editing = true;

  bool changed = false;
  if (ImGui::IsItemDeactivatedAfterEdit()) {
    st.editing = false;
    const char* t = st.buf;
    while (*t == ' ') ++t;
    if (*t == '\0') {
      if (cleared) *cleared = true;
    } else {
      float result = static_cast<float>(expr::eval(st.buf));
      if (!std::isnan(result)) {
        changed = (result != *v);
        *v = result;
      }
    }
    st.buf[0] = '\0';  // re-sync from *v on next idle frame
  }
  return changed;
}

static bool ExprInputInt(const char* id_str, int* v, bool* cleared = nullptr) {
  float fv = static_cast<float>(*v);
  bool  cl  = false;
  bool  ch  = ExprInputFloat(id_str, &fv, &cl);
  if (cl) { if (cleared) *cleared = true; return false; }
  if (ch) { *v = static_cast<int>(std::lround(static_cast<double>(fv))); return true; }
  return false;
}

// ============================================================
// UI helpers — standardised input style
// (label above, input box + inline unit label)
// ============================================================

static constexpr float kUnitColWidth = 44.0f;

static bool LabeledFloat(const char* label, const char* id, float* v,
                          const char* unit = nullptr,
                          const char* fmt  = "%.6g") {
  ImGui::TextUnformatted(label);
  float avail = ImGui::GetContentRegionAvail().x;
  float w = unit ? avail - kUnitColWidth : -1.0f;
  ImGui::SetNextItemWidth(w);
  bool ch = ExprInputFloat(id, v, nullptr, fmt);
  if (unit) { ImGui::SameLine(0, 6); ImGui::TextDisabled("%s", unit); }
  return ch;
}

static bool LabeledInt(const char* label, const char* id, int* v,
                        const char* unit  = nullptr) {
  ImGui::TextUnformatted(label);
  float avail = ImGui::GetContentRegionAvail().x;
  float w = unit ? avail - kUnitColWidth : -1.0f;
  ImGui::SetNextItemWidth(w);
  bool ch = ExprInputInt(id, v, nullptr);
  if (unit) { ImGui::SameLine(0, 6); ImGui::TextDisabled("%s", unit); }
  return ch;
}

// ============================================================
// Physics model definitions
// ============================================================

// Transient convection-diffusion-reaction with spatially-varying coefficients.
class ConvDiffModel : public mphys::Model {
 public:
  ConvDiffModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv,
                mphys::Field D_face, mphys::Field u_face,
                std::vector<double> k_cell,
                mphys::BoundaryCondition left_bc, mphys::BoundaryCondition right_bc)
      : Model(mesh, sv),
        D_face_(std::move(D_face)),
        u_face_(std::move(u_face)),
        k_cell_(std::move(k_cell)) {
    c_ = AddField("c", 0.0);
    SetBcs(c_, {left_bc, right_bc});
  }

  void Residual(double,
                const std::vector<mphys::Field>& y,
                const std::vector<mphys::Field>& ydot,
                const std::vector<double>&, const std::vector<double>&,
                std::vector<mphys::Field>& rr, std::vector<double>&) override {
    rr[c_] = mphys::fvm::Ddt(ydot[c_])
           + mphys::fvm::Convection(y[c_], u_face_, mesh_, bcs_[c_])
           - mphys::fvm::Laplacian(y[c_], D_face_, mesh_, bcs_[c_]);
    for (int i = 0; i < mesh_.n_cells; ++i)
      rr[c_][i] += y[c_][i] * k_cell_[i];
  }

 private:
  int c_ = 0;
  mphys::Field D_face_, u_face_;
  std::vector<double> k_cell_;
};

// Steady diffusion with spatially-varying face diffusivity.
class SteadyDiffModel : public mphys::Model {
 public:
  SteadyDiffModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv,
                  mphys::Field D_face,
                  mphys::BoundaryCondition left_bc, mphys::BoundaryCondition right_bc)
      : Model(mesh, sv), D_face_(std::move(D_face)) {
    double init = 0.5 * (left_bc.value + right_bc.value);
    c_ = AddField("c", init);
    SetBcs(c_, {left_bc, right_bc});
  }

  void Residual(double,
                const std::vector<mphys::Field>& y,
                const std::vector<mphys::Field>&,
                const std::vector<double>&, const std::vector<double>&,
                std::vector<mphys::Field>& rr, std::vector<double>&) override {
    rr[c_] = -mphys::fvm::Laplacian(y[c_], D_face_, mesh_, bcs_[c_]);
  }

 private:
  int c_ = 0;
  mphys::Field D_face_;
};

// 1D packed bed with Darcy's law + ideal-gas mass balance.
// Two coupled equations kept explicit:
//   (1) Darcy's Law  (algebraic): u + (κ/μ)·∂P/∂x = 0
//   (2) Mass balance (transient): ∂P/∂t + ∂(P·u)/∂x = 0
class DarcyPackedBedModel : public mphys::Model {
 public:
  DarcyPackedBedModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv,
                      std::vector<double> kappa,  // per-cell permeability [m²]
                      std::vector<double> mu,     // per-cell viscosity    [Pa·s]
                      const std::vector<double>& P_init,  // per-cell initial pressure
                      mphys::BoundaryCondition P_lbc, mphys::BoundaryCondition P_rbc,
                      mphys::BoundaryCondition u_lbc, mphys::BoundaryCondition u_rbc)
      : Model(mesh, sv), kappa_(std::move(kappa)), mu_(std::move(mu)) {
    P_ = AddField("P", 0.0);
    u_ = AddField("u", 0.0);
    auto& Pf = fields()[P_];
    for (int i = 0; i < mesh_.n_cells; ++i) Pf[i] = P_init[i];
    sv_.MarkFieldAlgebraic("u");   // Darcy has no ∂u/∂t — u is algebraic
    SetBcs(P_, {P_lbc, P_rbc});
    SetBcs(u_, {u_lbc, u_rbc});
  }

  void Residual(double,
                const std::vector<mphys::Field>& y,
                const std::vector<mphys::Field>& ydot,
                const std::vector<double>&, const std::vector<double>&,
                std::vector<mphys::Field>& rr, std::vector<double>&) override {
    // ── Darcy's Law (algebraic) ──────────────────────────────────────────────
    //   u + (κ/μ) · ∂P/∂x = 0
    auto grad_P = mphys::fvm::Grad(y[P_], mesh_, bcs_[P_]);
    for (int i = 0; i < mesh_.n_cells; ++i) {
      double dPdx = 0.5 * (grad_P[i] + grad_P[i + 1]);  // cell-centre gradient
      rr[u_][i] = y[u_][i] + (kappa_[i] / mu_[i]) * dPdx;
    }
    // ── Ideal-gas mass balance (ρ = PM/RT, T constant → ρ ∝ P) ─────────────
    //   ∂P/∂t + ∂(P·u)/∂x = 0
    auto u_face = mphys::fvm::InterpolateToFaces(y[u_], mesh_, bcs_[u_]);
    rr[P_] = mphys::fvm::Ddt(ydot[P_])
           + mphys::fvm::Convection(y[P_], u_face, mesh_, bcs_[P_]);
  }

 private:
  int P_ = 0, u_ = 1;
  std::vector<double> kappa_, mu_;
};

// ============================================================
// Geometry and application state
// ============================================================

struct GeoNode   { float x = 0.0f; };

struct GeoDomain {
  int   n_cells = 20;
  // Convection-Diffusion-Reaction
  float D     = 1.0f;    // diffusivity   [m²/s]
  float u     = 0.0f;    // velocity      [m/s]
  float k     = 0.0f;    // reaction rate [1/s]
  // Darcy Packed Bed
  float kappa = 1e-10f;  // permeability  [m²]
  float mu    = 1.8e-5f; // dyn. viscosity [Pa·s]
};

struct NodeBc {
  int   type  = 0;    // 0 = Dirichlet, 1 = Neumann
  float value = 0.0f;
};

struct PackedBedBc {
  int   type  = 0;    // 0 = Pressure, 1 = Velocity
  float value = 0.0f;
};

struct Geometry1D {
  std::vector<GeoNode>   nodes;
  std::vector<GeoDomain> domains;
  bool built          = false;
  int  selected_node   = -1;
  int  selected_domain = -1;
};

enum class NavNode { Geometry, Physics, Mesh, Study, Results };
enum class PhysicsModule { ConvectionDiffusion, SteadyDiffusion, DarcyPackedBed,
                           SingleParticle, SingleParticleElectrolyte };

// User-facing Single Particle Model inputs (floats for ImGui widgets).
// Defaults are representative LG M50 (Chen 2020) values.
struct SpmInputs {
  // Geometry
  float R_n = 5.86e-6f, R_p = 5.22e-6f;   // particle radii          [m]
  float L_n = 85.2e-6f, L_p = 75.6e-6f;   // electrode thicknesses   [m]
  float A   = 0.1027f;                     // electrode area          [m²]
  float eps_n = 0.75f, eps_p = 0.665f;     // active-material fraction [-]
  // Transport
  float D_n = 3.3e-14f, D_p = 4.0e-15f;    // solid diffusivities     [m²/s]
  // Concentrations
  float cn_max = 33133.0f, cp_max = 63104.0f;  // max concentrations  [mol/m³]
  float x0 = 0.84f, y0 = 0.27f;            // initial stoichiometries [-]
  float c_e = 1000.0f;                      // electrolyte conc.      [mol/m³]
  // Kinetics
  float k_n = 6.48e-7f, k_p = 3.42e-6f;    // reaction rate constants
  // Operating
  float I = 5.0f, T = 298.15f;             // current [A], temperature [K]
  int   n_cells = 30;                       // mesh points per particle
  // Selected database materials (indices into the domain-filtered catalogues;
  // see MaterialId helpers).  The electrode choice supplies the OCP curve and
  // pre-fills the numeric fields above; the electrolyte choice the salt conc.
  int   neg_material   = 0;                 // negative electrode
  int   pos_material   = 0;                 // positive electrode
  int   elyte_material = 0;                 // electrolyte
};

// --- Database material <-> GUI helpers ------------------------------------
// The combos store an index into the (deterministic) domain-filtered id lists;
// these map an index back to the strongly-typed identifier.
static mphys::materials::ElectrodeId NegElectrodeId(int idx) {
  const auto ids =
      mphys::materials::Database::ElectrodeIds(mphys::materials::Domain::kNegativeElectrode);
  if (idx < 0 || idx >= static_cast<int>(ids.size())) idx = 0;
  return ids[idx];
}
static mphys::materials::ElectrodeId PosElectrodeId(int idx) {
  const auto ids =
      mphys::materials::Database::ElectrodeIds(mphys::materials::Domain::kPositiveElectrode);
  if (idx < 0 || idx >= static_cast<int>(ids.size())) idx = 0;
  return ids[idx];
}
static mphys::materials::ElectrolyteId ElectrolyteIdAt(int idx) {
  const auto ids = mphys::materials::Database::ElectrolyteIds();
  if (idx < 0 || idx >= static_cast<int>(ids.size())) idx = 0;
  return ids[idx];
}

// Load a material's scalar properties into the editable input fields.
static void ApplyNegElectrode(SpmInputs& in) {
  const auto& m = mphys::materials::Database::Electrode(NegElectrodeId(in.neg_material));
  in.cn_max = static_cast<float>(m.c_max);
  in.D_n    = static_cast<float>(m.diffusivity);
  in.k_n    = static_cast<float>(m.reaction_rate);
}
static void ApplyPosElectrode(SpmInputs& in) {
  const auto& m = mphys::materials::Database::Electrode(PosElectrodeId(in.pos_material));
  in.cp_max = static_cast<float>(m.c_max);
  in.D_p    = static_cast<float>(m.diffusivity);
  in.k_p    = static_cast<float>(m.reaction_rate);
}

// Full-width dropdown over a list of display names. Returns true if the
// selection changed (and writes the new index into *idx).
static bool MaterialCombo(const char* combo_id,
                          const std::vector<std::string>& names, int* idx) {
  if (*idx < 0 || *idx >= static_cast<int>(names.size())) *idx = 0;
  std::vector<const char*> items;
  items.reserve(names.size());
  for (const auto& n : names) items.push_back(n.c_str());
  ImGui::SetNextItemWidth(-1.0f);
  return ImGui::Combo(combo_id, idx, items.data(),
                      static_cast<int>(items.size()));
}

static mphys::models::SpmParameters ToSpmParameters(const SpmInputs& in) {
  mphys::models::SpmParameters p;
  p.R_n = in.R_n;     p.R_p = in.R_p;
  p.L_n = in.L_n;     p.L_p = in.L_p;
  p.A   = in.A;
  p.eps_n = in.eps_n; p.eps_p = in.eps_p;
  p.D_n = in.D_n;     p.D_p = in.D_p;
  p.cn_max = in.cn_max; p.cp_max = in.cp_max;
  p.cn0 = in.x0 * in.cn_max;
  p.cp0 = in.y0 * in.cp_max;
  p.c_e = in.c_e;
  p.k_n = in.k_n;     p.k_p = in.k_p;
  p.I = in.I;         p.T = in.T;
  // Open-circuit potentials come from the selected electrode materials — they
  // are curves, not editable scalars, so the database is their only source
  // (this is what makes e.g. an LFP positive electrode behave correctly).
  p.Un = mphys::materials::Database::Electrode(NegElectrodeId(in.neg_material)).ocp;
  p.Up = mphys::materials::Database::Electrode(PosElectrodeId(in.pos_material)).ocp;
  return p;
}

// User-facing SPMe inputs: the SPM particle/kinetics core plus the macroscopic
// electrolyte (separator + porous-electrode electrolyte transport) parameters.
struct SpmeInputs {
  SpmInputs core;                         // particle, kinetics, operating
  // Separator
  float L_s = 12e-6f;                      // separator thickness     [m]
  // Electrolyte volume fractions (porosity) per region
  float eps_e_n = 0.25f, eps_e_s = 0.47f, eps_e_p = 0.335f;
  // Electrolyte transport / thermodynamics
  float D_e     = 1.769e-10f;              // salt diffusivity        [m²/s]
  float kappa_e = 1.0f;                    // ionic conductivity      [S/m]
  float t_plus  = 0.2594f;                 // transference number     [-]
  float brugg   = 1.5f;                     // Bruggeman exponent      [-]
  // Solid conductivities
  float sigma_n = 215.0f, sigma_p = 0.18f; // solid conductivities    [S/m]
  // Initial electrolyte concentration
  float ce0 = 1000.0f;                     // initial conc.           [mol/m³]
  // Mesh: cells per region (the particle mesh uses the sum)
  int n_n = 20, n_s = 12, n_p = 20;
};

static mphys::models::SpmeParameters ToSpmeParameters(const SpmeInputs& in) {
  mphys::models::SpmeParameters p;
  p.core    = ToSpmParameters(in.core);
  p.L_s     = in.L_s;
  p.eps_e_n = in.eps_e_n; p.eps_e_s = in.eps_e_s; p.eps_e_p = in.eps_e_p;
  p.D_e     = in.D_e;     p.kappa_e = in.kappa_e;
  p.t_plus  = in.t_plus;  p.brugg   = in.brugg;
  p.sigma_n = in.sigma_n; p.sigma_p = in.sigma_p;
  p.ce0     = in.ce0;
  return p;
}

struct AppState {
  NavNode       nav     = NavNode::Geometry;
  PhysicsModule physics = PhysicsModule::ConvectionDiffusion;

  mphys::CoordSystem coord_system = mphys::CoordSystem::kCartesian;

  Geometry1D geo = {
    {{0.0f}, {1.0f}},
    {{20, 1.0f, 0.0f, 0.0f}},
    false, -1, -1
  };

  int geo_input_mode = 0;                        // 0=coordinates, 1=lengths
  std::vector<float> geo_lengths = {1.0f};
  int geo_unit = 0;                              // 0=m  1=cm  2=mm  3=µm

  NodeBc left_bc  = {0, 1.0f};
  NodeBc right_bc = {1, 0.0f};
  // Darcy packed bed BCs (pressure or velocity at each terminal node)
  PackedBedBc pb_left_bc  = {0, 2.0f};   // Pressure = 2 (high)
  PackedBedBc pb_right_bc = {0, 1.0f};   // Pressure = 1 (low)

  float t_end        = 50.0f;
  float dt_initial   = 1e-3f;
  float dt_max       = 1.0f;
  float dt_snapshot  = 1.0f;   // interval at which snapshots are recorded
  float rel_tol      = 1e-6f;
  float abs_tol      = 1e-8f;

  SpmInputs            spm;
  SpmeInputs           spme;
  std::vector<double>  spm_time;        // voltage-curve time axis [s]
  std::vector<double>  spm_voltage;     // terminal voltage        [V]
  std::vector<double>  spme_ce_x;       // electrolyte cell centres [m] (x-axis)

  mphys::SimResult     result;
  std::vector<double>  cell_centres;
  float                plot_time    = 0.0f;
  bool                 has_results  = false;
  std::string          status_msg;

  bool dark_theme = true;
};

// ============================================================
// Face-field computation helpers
// ============================================================

// Returns cell-domain index for every cell in the composite mesh.
static std::vector<int> CellDomainMap(const Geometry1D& geo) {
  std::vector<int> m;
  m.reserve(0);
  for (int d = 0; d < (int)geo.domains.size(); ++d)
    for (int i = 0; i < geo.domains[d].n_cells; ++i)
      m.push_back(d);
  return m;
}

// Harmonic-mean face diffusivity accounting for heterogeneous dx.
// D_face[i] for interior face i uses cells i-1 (left) and i (right).
static mphys::Field ComputeFaceD(const Geometry1D& geo, const mphys::Mesh1D& mesh) {
  auto dom = CellDomainMap(geo);
  int  n   = mesh.n_cells;
  mphys::Field D_face("D_face", n + 1);
  D_face[0] = geo.domains[dom[0]].D;
  D_face[n] = geo.domains[dom[n - 1]].D;
  for (int i = 1; i < n; ++i) {
    double D_L  = geo.domains[dom[i - 1]].D;
    double D_R  = geo.domains[dom[i]].D;
    double dx_L = mesh.dx[i - 1];
    double dx_R = mesh.dx[i];
    // Resistance-weighted harmonic mean: D = (dxL+dxR) / (dxL/DL + dxR/DR)
    D_face[i] = (dx_L + dx_R) / (dx_L / D_L + dx_R / D_R);
  }
  return D_face;
}

// Distance-weighted arithmetic-mean face velocity.
static mphys::Field ComputeFaceU(const Geometry1D& geo, const mphys::Mesh1D& mesh) {
  auto dom = CellDomainMap(geo);
  int  n   = mesh.n_cells;
  mphys::Field u_face("u_face", n + 1);
  u_face[0] = geo.domains[dom[0]].u;
  u_face[n] = geo.domains[dom[n - 1]].u;
  for (int i = 1; i < n; ++i) {
    double u_L  = geo.domains[dom[i - 1]].u;
    double u_R  = geo.domains[dom[i]].u;
    double dx_L = mesh.dx[i - 1];
    double dx_R = mesh.dx[i];
    u_face[i] = (u_R * dx_L + u_L * dx_R) / (dx_L + dx_R);
  }
  return u_face;
}

// Cell-centred reaction rate.
static std::vector<double> ComputeCellK(const Geometry1D& geo) {
  std::vector<double> k;
  for (const auto& d : geo.domains)
    for (int i = 0; i < d.n_cells; ++i)
      k.push_back(d.k);
  return k;
}

static std::vector<double> ComputeCellKappa(const Geometry1D& geo) {
  std::vector<double> v;
  for (const auto& d : geo.domains)
    for (int i = 0; i < d.n_cells; ++i)
      v.push_back(d.kappa);
  return v;
}

static std::vector<double> ComputeCellMu(const Geometry1D& geo) {
  std::vector<double> v;
  for (const auto& d : geo.domains)
    for (int i = 0; i < d.n_cells; ++i)
      v.push_back(d.mu);
  return v;
}

// ============================================================
// Composite mesh builder
// ============================================================

static mphys::Mesh1D BuildCompositeMesh(const Geometry1D& geo, mphys::CoordSystem cs) {
  mphys::Mesh1D mesh;
  for (int d = 0; d < (int)geo.domains.size(); ++d) {
    double x0 = geo.nodes[d].x;
    double x1 = geo.nodes[d + 1].x;
    int    nc  = geo.domains[d].n_cells;
    auto   seg = mphys::MakeUniformMesh1D(x0, x1, nc, cs);
    for (int i = (d == 0 ? 0 : 1); i <= nc; ++i)
      mesh.face_positions.push_back(seg.face_positions[i]);
    for (auto v : seg.cell_centres) mesh.cell_centres.push_back(v);
    for (auto v : seg.dx)           mesh.dx.push_back(v);
  }
  mesh.n_cells = (int)mesh.cell_centres.size();
  mesh.df.resize(mesh.n_cells + 1);
  mesh.df[0] = 2.0 * (mesh.cell_centres[0] - mesh.face_positions[0]);
  for (int i = 1; i < mesh.n_cells; ++i)
    mesh.df[i] = mesh.cell_centres[i] - mesh.cell_centres[i - 1];
  mesh.df[mesh.n_cells] = 2.0 * (mesh.face_positions.back() - mesh.cell_centres.back());
  return mesh;
}

// ============================================================
// Simulation runner
// ============================================================

// Single Particle Model: two normalised spherical particles + voltage output.
static void RunSpm(AppState& s) {
  try {
    auto p = ToSpmParameters(s.spm);
    int nc = std::max(2, s.spm.n_cells);

    auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, nc,
                                         mphys::CoordSystem::kSpherical);
    s.cell_centres = mesh.cell_centres;   // normalised radius r/R ∈ [0,1]

    mphys::StateVector sv(mesh.n_cells);
    mphys::models::SpmModel model(mesh, sv, p);

    mphys::SolverOptions opts;
    opts.tolerance.relative = static_cast<double>(s.rel_tol);
    opts.tolerance.absolute = static_cast<double>(s.abs_tol);
    opts.initial_time_step  = static_cast<double>(s.dt_initial);
    opts.maximum_time_step  = static_cast<double>(s.dt_max);

    double dt_snap = std::max(static_cast<double>(s.dt_snapshot),
                              static_cast<double>(s.dt_max));

    mphys::SunContext sunctx;
    mphys::TransientSolver solver(model, opts, sunctx);
    double next_snap = 0.0;
    std::string warn = solver.Solve(0.0, static_cast<double>(s.t_end),
        [&](double t, const std::vector<mphys::Field>& f,
            const std::vector<double>& a) {
          if (t >= next_snap - 1e-12) {
            s.result.Record(t, f, a);
            s.spm_time.push_back(t);
            s.spm_voltage.push_back(a[model.voltage_index()]);
            next_snap += dt_snap;
          }
        });

    if (!s.result.snapshots.empty()) {
      s.plot_time   = static_cast<float>(s.result.snapshots.back().t);
      s.has_results = true;
      s.nav         = NavNode::Results;
    }
    s.status_msg = warn.empty()
        ? "Done — " + std::to_string(s.result.snapshots.size()) + " snapshots"
        : "Warning: " + warn;
  } catch (const std::exception& e) {
    s.status_msg = std::string("Error: ") + e.what();
  }
}

// SPMe: two spherical particles + a macroscopic electrolyte field, plus the
// terminal-voltage output.  Particles and electrolyte share a cell count but
// live on different meshes (spherical r∈[0,1] vs Cartesian x∈[0,L]).
static void RunSpme(AppState& s) {
  try {
    auto p = ToSpmeParameters(s.spme);
    int nn = std::max(2, s.spme.n_n);
    int ns = std::max(1, s.spme.n_s);
    int np = std::max(2, s.spme.n_p);

    auto em = mphys::models::MakeSpmeElectrolyteMesh(p, nn, ns, np);
    auto particle_mesh = mphys::MakeUniformMesh1D(
        0.0, 1.0, em.mesh.n_cells, mphys::CoordSystem::kSpherical);

    s.cell_centres = particle_mesh.cell_centres;  // normalised radius r/R for particles
    s.spme_ce_x    = em.mesh.cell_centres;         // electrolyte position x [m]

    mphys::StateVector sv(em.mesh.n_cells);
    mphys::models::SpmeModel model(particle_mesh, em, sv, p);

    mphys::SolverOptions opts;
    opts.tolerance.relative = static_cast<double>(s.rel_tol);
    opts.tolerance.absolute = static_cast<double>(s.abs_tol);
    opts.initial_time_step  = static_cast<double>(s.dt_initial);
    opts.maximum_time_step  = static_cast<double>(s.dt_max);

    double dt_snap = std::max(static_cast<double>(s.dt_snapshot),
                              static_cast<double>(s.dt_max));

    mphys::SunContext sunctx;
    mphys::TransientSolver solver(model, opts, sunctx);
    double next_snap = 0.0;
    std::string warn = solver.Solve(0.0, static_cast<double>(s.t_end),
        [&](double t, const std::vector<mphys::Field>& f,
            const std::vector<double>& a) {
          if (t >= next_snap - 1e-12) {
            s.result.Record(t, f, a);
            s.spm_time.push_back(t);
            s.spm_voltage.push_back(a[model.voltage_index()]);
            next_snap += dt_snap;
          }
        });

    if (!s.result.snapshots.empty()) {
      s.plot_time   = static_cast<float>(s.result.snapshots.back().t);
      s.has_results = true;
      s.nav         = NavNode::Results;
    }
    s.status_msg = warn.empty()
        ? "Done — " + std::to_string(s.result.snapshots.size()) + " snapshots"
        : "Warning: " + warn;
  } catch (const std::exception& e) {
    s.status_msg = std::string("Error: ") + e.what();
  }
}

static void RunSimulation(AppState& s) {
  s.has_results = false;
  s.result.snapshots.clear();
  s.cell_centres.clear();
  s.spme_ce_x.clear();
  s.spm_time.clear();
  s.spm_voltage.clear();
  s.status_msg = "Running...";

  if (s.physics == PhysicsModule::SingleParticle) { RunSpm(s); return; }
  if (s.physics == PhysicsModule::SingleParticleElectrolyte) { RunSpme(s); return; }

  try {
    if (!s.geo.built || s.geo.domains.empty())
      throw std::runtime_error("Build the geometry first (Geometry -> Build)");
    for (int d = 0; d < (int)s.geo.domains.size(); ++d)
      if (s.geo.nodes[d + 1].x <= s.geo.nodes[d].x)
        throw std::runtime_error("Node positions must be strictly increasing");

    auto mesh  = BuildCompositeMesh(s.geo, s.coord_system);
    s.cell_centres = mesh.cell_centres;

    auto D_face = ComputeFaceD(s.geo, mesh);
    auto u_face = ComputeFaceU(s.geo, mesh);
    auto k_cell = ComputeCellK(s.geo);

    auto make_bc = [](const NodeBc& nb) -> mphys::BoundaryCondition {
      return nb.type == 0 ? mphys::DirichletBc(nb.value) : mphys::NeumannBc(nb.value);
    };
    auto lbc = make_bc(s.left_bc);
    auto rbc = make_bc(s.right_bc);

    mphys::SolverOptions opts;
    opts.tolerance.relative = static_cast<double>(s.rel_tol);
    opts.tolerance.absolute = static_cast<double>(s.abs_tol);
    opts.initial_time_step  = static_cast<double>(s.dt_initial);
    opts.maximum_time_step  = static_cast<double>(s.dt_max);

    mphys::SunContext  sunctx;
    mphys::StateVector sv(mesh.n_cells);

    double dt_snap = std::max(static_cast<double>(s.dt_snapshot),
                               static_cast<double>(s.dt_max));
    std::string solver_warning;
    if (s.physics == PhysicsModule::DarcyPackedBed) {
      auto kappa_cell = ComputeCellKappa(s.geo);
      auto mu_cell    = ComputeCellMu(s.geo);

      double kL = s.geo.domains.front().kappa, muL = s.geo.domains.front().mu;
      double kR = s.geo.domains.back().kappa,  muR = s.geo.domains.back().mu;

      // Pressure BC → Dirichlet on P; Velocity BC → Neumann on P (from Darcy)
      auto make_P_bc = [](const PackedBedBc& bc, double kp, double mp) {
        return bc.type == 0 ? mphys::DirichletBc(bc.value)
                            : mphys::NeumannBc(-mp * bc.value / kp);
      };
      // Velocity BC → Dirichlet on u; Pressure BC → Neumann on u (free)
      auto make_u_bc = [](const PackedBedBc& bc) {
        return bc.type == 1 ? mphys::DirichletBc(bc.value)
                            : mphys::NeumannBc(0.0);
      };

      auto P_lbc = make_P_bc(s.pb_left_bc,  kL, muL);
      auto P_rbc = make_P_bc(s.pb_right_bc, kR, muR);
      auto u_lbc = make_u_bc(s.pb_left_bc);
      auto u_rbc = make_u_bc(s.pb_right_bc);

      // Initialise P with a linear profile so IDACalcIC has a consistent starting point.
      double P_scalar = (s.pb_left_bc.type == 0) ? s.pb_left_bc.value
                      : (s.pb_right_bc.type == 0) ? s.pb_right_bc.value
                      : 101325.0;
      double P_L_val  = (s.pb_left_bc.type  == 0) ? (double)s.pb_left_bc.value  : P_scalar;
      double P_R_val  = (s.pb_right_bc.type == 0) ? (double)s.pb_right_bc.value : P_scalar;
      std::vector<double> P_init_vec(mesh.n_cells);
      {
        double x_L = mesh.cell_centres.front(), x_R = mesh.cell_centres.back();
        double span = (x_R > x_L) ? (x_R - x_L) : 1.0;
        for (int i = 0; i < mesh.n_cells; ++i)
          P_init_vec[i] = P_L_val + (P_R_val - P_L_val) *
                          (mesh.cell_centres[i] - x_L) / span;
      }

      DarcyPackedBedModel model(mesh, sv, kappa_cell, mu_cell, P_init_vec,
                                P_lbc, P_rbc, u_lbc, u_rbc);
      mphys::TransientSolver solver(model, opts, sunctx);
      double next_snap = 0.0;
      solver_warning = solver.Solve(0.0, static_cast<double>(s.t_end),
          [&](double t, const std::vector<mphys::Field>& f,
              const std::vector<double>& a) {
            if (t >= next_snap - 1e-12) {
              s.result.Record(t, f, a);
              next_snap += dt_snap;
            }
          });
    } else if (s.physics == PhysicsModule::ConvectionDiffusion) {
      ConvDiffModel model(mesh, sv, D_face, u_face, k_cell, lbc, rbc);
      mphys::TransientSolver solver(model, opts, sunctx);
      double next_snap = 0.0;
      solver_warning = solver.Solve(0.0, static_cast<double>(s.t_end),
          [&](double t, const std::vector<mphys::Field>& f,
              const std::vector<double>& a) {
            if (t >= next_snap - 1e-12) {
              s.result.Record(t, f, a);
              next_snap += dt_snap;
            }
          });
    } else {
      SteadyDiffModel model(mesh, sv, D_face, lbc, rbc);
      mphys::SteadySolver solver(model, opts, sunctx);
      solver.Solve();
      s.result.Record(0.0, model.fields(), model.algebraics());
    }

    if (!s.result.snapshots.empty()) {
      s.plot_time   = static_cast<float>(s.result.snapshots.back().t);
      s.has_results = true;
      s.nav         = NavNode::Results;
    }
    if (!solver_warning.empty())
      s.status_msg = "Warning: " + solver_warning +
                     " (" + std::to_string(s.result.snapshots.size()) + " snapshots)";
    else
      s.status_msg = "Done — " + std::to_string(s.result.snapshots.size()) + " snapshots";

  } catch (const std::exception& e) {
    s.status_msg = std::string("Error: ") + e.what();
  }
}

// ============================================================
// Cereal serialisation — non-intrusive, JSON-friendly
// ============================================================

template<class Ar> void serialize(Ar& ar, GeoNode& v) {
  ar(cereal::make_nvp("x", v.x));
}

template<class Ar> void serialize(Ar& ar, GeoDomain& v) {
  ar(cereal::make_nvp("n_cells", v.n_cells),
     cereal::make_nvp("D",       v.D),
     cereal::make_nvp("u",       v.u),
     cereal::make_nvp("k",       v.k),
     cereal::make_nvp("kappa",   v.kappa),
     cereal::make_nvp("mu",      v.mu));
}

template<class Ar> void serialize(Ar& ar, NodeBc& v) {
  ar(cereal::make_nvp("type",  v.type),
     cereal::make_nvp("value", v.value));
}

template<class Ar> void serialize(Ar& ar, PackedBedBc& v) {
  ar(cereal::make_nvp("type",  v.type),
     cereal::make_nvp("value", v.value));
}

template<class Ar> void serialize(Ar& ar, Geometry1D& v) {
  ar(cereal::make_nvp("nodes",   v.nodes),
     cereal::make_nvp("domains", v.domains),
     cereal::make_nvp("built",   v.built));
}

template<class Ar> void serialize(Ar& ar, SpmInputs& v) {
  ar(cereal::make_nvp("R_n", v.R_n), cereal::make_nvp("R_p", v.R_p),
     cereal::make_nvp("L_n", v.L_n), cereal::make_nvp("L_p", v.L_p),
     cereal::make_nvp("A", v.A),
     cereal::make_nvp("eps_n", v.eps_n), cereal::make_nvp("eps_p", v.eps_p),
     cereal::make_nvp("D_n", v.D_n), cereal::make_nvp("D_p", v.D_p),
     cereal::make_nvp("cn_max", v.cn_max), cereal::make_nvp("cp_max", v.cp_max),
     cereal::make_nvp("x0", v.x0), cereal::make_nvp("y0", v.y0),
     cereal::make_nvp("c_e", v.c_e),
     cereal::make_nvp("k_n", v.k_n), cereal::make_nvp("k_p", v.k_p),
     cereal::make_nvp("I", v.I), cereal::make_nvp("T", v.T),
     cereal::make_nvp("n_cells", v.n_cells),
     cereal::make_nvp("neg_material", v.neg_material),
     cereal::make_nvp("pos_material", v.pos_material),
     cereal::make_nvp("elyte_material", v.elyte_material));
}

template<class Ar> void serialize(Ar& ar, SpmeInputs& v) {
  ar(cereal::make_nvp("core", v.core),
     cereal::make_nvp("L_s", v.L_s),
     cereal::make_nvp("eps_e_n", v.eps_e_n), cereal::make_nvp("eps_e_s", v.eps_e_s),
     cereal::make_nvp("eps_e_p", v.eps_e_p),
     cereal::make_nvp("D_e", v.D_e), cereal::make_nvp("kappa_e", v.kappa_e),
     cereal::make_nvp("t_plus", v.t_plus), cereal::make_nvp("brugg", v.brugg),
     cereal::make_nvp("sigma_n", v.sigma_n), cereal::make_nvp("sigma_p", v.sigma_p),
     cereal::make_nvp("ce0", v.ce0),
     cereal::make_nvp("n_n", v.n_n), cereal::make_nvp("n_s", v.n_s),
     cereal::make_nvp("n_p", v.n_p));
}

// Split save/load for AppState so enums round-trip cleanly as ints.
template<class Ar> void save(Ar& ar, const AppState& v) {
  ar(cereal::make_nvp("physics",        static_cast<int>(v.physics)),
     cereal::make_nvp("coord_system",   static_cast<int>(v.coord_system)),
     cereal::make_nvp("geo",            v.geo),
     cereal::make_nvp("geo_input_mode", v.geo_input_mode),
     cereal::make_nvp("geo_lengths",    v.geo_lengths),
     cereal::make_nvp("geo_unit",       v.geo_unit),
     cereal::make_nvp("left_bc",        v.left_bc),
     cereal::make_nvp("right_bc",       v.right_bc),
     cereal::make_nvp("pb_left_bc",     v.pb_left_bc),
     cereal::make_nvp("pb_right_bc",    v.pb_right_bc),
     cereal::make_nvp("spm",            v.spm),
     cereal::make_nvp("spme",           v.spme),
     cereal::make_nvp("t_end",          v.t_end),
     cereal::make_nvp("dt_initial",     v.dt_initial),
     cereal::make_nvp("dt_max",         v.dt_max),
     cereal::make_nvp("dt_snapshot",    v.dt_snapshot),
     cereal::make_nvp("rel_tol",        v.rel_tol),
     cereal::make_nvp("abs_tol",        v.abs_tol),
     cereal::make_nvp("dark_theme",     v.dark_theme));
}

template<class Ar> void load(Ar& ar, AppState& v) {
  int physics_i = 0, coord_i = 0;
  ar(cereal::make_nvp("physics",        physics_i),
     cereal::make_nvp("coord_system",   coord_i),
     cereal::make_nvp("geo",            v.geo),
     cereal::make_nvp("geo_input_mode", v.geo_input_mode),
     cereal::make_nvp("geo_lengths",    v.geo_lengths),
     cereal::make_nvp("geo_unit",       v.geo_unit),
     cereal::make_nvp("left_bc",        v.left_bc),
     cereal::make_nvp("right_bc",       v.right_bc),
     cereal::make_nvp("pb_left_bc",     v.pb_left_bc),
     cereal::make_nvp("pb_right_bc",    v.pb_right_bc),
     cereal::make_nvp("spm",            v.spm),
     cereal::make_nvp("spme",           v.spme),
     cereal::make_nvp("t_end",          v.t_end),
     cereal::make_nvp("dt_initial",     v.dt_initial),
     cereal::make_nvp("dt_max",         v.dt_max),
     cereal::make_nvp("dt_snapshot",    v.dt_snapshot),
     cereal::make_nvp("rel_tol",        v.rel_tol),
     cereal::make_nvp("abs_tol",        v.abs_tol),
     cereal::make_nvp("dark_theme",     v.dark_theme));
  v.physics      = static_cast<PhysicsModule>(physics_i);
  v.coord_system = static_cast<mphys::CoordSystem>(coord_i);
}

static void SaveState(const AppState& s, const std::string& path) {
  std::ofstream os(path);
  if (!os) throw std::runtime_error("Cannot open for writing: " + path);
  cereal::JSONOutputArchive ar(os);
  ar(cereal::make_nvp("mphys_state", const_cast<AppState&>(s)));
}

static void LoadState(AppState& s, const std::string& path) {
  std::ifstream is(path);
  if (!is) throw std::runtime_error("Cannot open file: " + path);
  {
    cereal::JSONInputArchive ar(is);
    ar(cereal::make_nvp("mphys_state", s));
  }
  // Clear runtime-only state
  s.nav         = NavNode::Geometry;
  s.has_results = false;
  s.result.snapshots.clear();
  s.cell_centres.clear();
  s.spme_ce_x.clear();
  s.plot_time   = 0.0f;
  s.status_msg.clear();
  s.geo.selected_node = s.geo.selected_domain = -1;
}

// ============================================================
// Geometry helpers
// ============================================================

static void ApplyLengthsToNodes(AppState& s) {
  s.geo.nodes.clear();
  s.geo.nodes.push_back({0.0f});
  float x = 0.0f;
  for (float L : s.geo_lengths) { x += L; s.geo.nodes.push_back({x}); }
  while ((int)s.geo.domains.size() < (int)s.geo_lengths.size())
    s.geo.domains.push_back({});
  s.geo.domains.resize(s.geo_lengths.size());
}

static void ApplyNodesToLengths(AppState& s) {
  s.geo_lengths.clear();
  for (int d = 0; d < (int)s.geo.domains.size(); ++d)
    s.geo_lengths.push_back(s.geo.nodes[d + 1].x - s.geo.nodes[d].x);
}

// ============================================================
// GUI panels
// ============================================================

static void ShowGeometryPanel(AppState& s) {
  static constexpr const char* kUnits[] = {"m", "cm", "mm", "\xc2\xb5m"};
  Geometry1D& geo = s.geo;

  if (s.physics == PhysicsModule::SingleParticle) {
    ImGui::SeparatorText("Single Particle Model");
    ImGui::Spacing();
    ImGui::TextWrapped(
        "The SPM builds its own geometry: two spherical particles on a "
        "normalised radius r/R in [0,1]. No domain setup is required.");
    ImGui::Spacing();
    ImGui::BulletText("Set particle sizes & materials in the Physics node.");
    ImGui::BulletText("Set radial resolution in the Mesh node.");
    ImGui::BulletText("Set current & end time in the Study node, then Run.");
    return;
  }

  if (s.physics == PhysicsModule::SingleParticleElectrolyte) {
    ImGui::SeparatorText("Single Particle Model w/ Electrolyte");
    ImGui::Spacing();
    ImGui::TextWrapped(
        "The SPMe builds its own geometry: two spherical particles (r/R in "
        "[0,1]) plus a macroscopic electrolyte field on the negative | "
        "separator | positive sandwich. No domain setup is required.");
    ImGui::Spacing();
    ImGui::BulletText("Set particle, electrolyte & separator data in Physics.");
    ImGui::BulletText("Set cells per region in the Mesh node.");
    ImGui::BulletText("Set current & end time in the Study node, then Run.");
    return;
  }

  ImGui::SeparatorText("Coordinate System");
  ImGui::Spacing();
  int cs_sel = (s.coord_system == mphys::CoordSystem::kSpherical) ? 1 : 0;
  if (ImGui::RadioButton("Cartesian", &cs_sel, 0))
    s.coord_system = mphys::CoordSystem::kCartesian;
  ImGui::SameLine();
  if (ImGui::RadioButton("Spherical", &cs_sel, 1))
    s.coord_system = mphys::CoordSystem::kSpherical;
  ImGui::SameLine(0, 20.0f);
  ImGui::TextUnformatted("Unit:");
  ImGui::SameLine();
  ImGui::SetNextItemWidth(60.0f);
  ImGui::Combo("##unit", &s.geo_unit, kUnits, 4);

  ImGui::Spacing();
  ImGui::SeparatorText("Domain Definition");
  ImGui::Spacing();

  int old_mode = s.geo_input_mode;
  ImGui::TextUnformatted("Specify by:");
  ImGui::SameLine();
  ImGui::RadioButton("Coordinates##m", &s.geo_input_mode, 0);
  ImGui::SameLine();
  ImGui::RadioButton("Lengths##m",     &s.geo_input_mode, 1);

  if (s.geo_input_mode != old_mode) {
    if (s.geo_input_mode == 1) ApplyNodesToLengths(s);
    else                        ApplyLengthsToNodes(s);
    geo.built = false;
  }

  ImGui::Spacing();

  // Lock origin to 0
  if (!geo.nodes.empty()) geo.nodes[0].x = 0.0f;

  static constexpr ImGuiTableFlags kTblFlags =
      ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH |
      ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoHostExtendX;

  if (s.geo_input_mode == 0) {
    char hint[64]; snprintf(hint, sizeof(hint),
        "x positions  [%s]  (clear a row to remove it)", kUnits[s.geo_unit]);
    ImGui::TextDisabled("%s", hint);
    ImGui::Spacing();

    if (ImGui::BeginTable("coord_tbl", 1, kTblFlags)) {
      int delete_idx = -1;

      for (int i = 0; i < (int)geo.nodes.size(); ++i) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::SetNextItemWidth(-1.0f);
        char id[32]; snprintf(id, sizeof(id), "##cx%d", i);

        if (i == 0) {
          ImGui::BeginDisabled();
          float zero = 0.0f;
          ImGui::InputFloat(id, &zero, 0, 0, "%.6g");
          ImGui::EndDisabled();
        } else {
          bool cleared = false;
          if (ExprInputFloat(id, &geo.nodes[i].x, &cleared)) geo.built = false;
          if (cleared && (int)geo.nodes.size() > 2) delete_idx = i;
        }
      }

      if (delete_idx >= 0) {
        int di = delete_idx - 1;
        geo.nodes.erase(geo.nodes.begin() + delete_idx);
        if (di < (int)geo.domains.size())
          geo.domains.erase(geo.domains.begin() + di);
        geo.built = false;
        geo.selected_node = geo.selected_domain = -1;
        ApplyNodesToLengths(s);
      }

      {
        static char coord_draft_buf[32] = "";
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::SetNextItemWidth(-1.0f);
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyle().Colors[ImGuiCol_TextDisabled]);
        ImGui::InputText("##cdraft", coord_draft_buf, sizeof(coord_draft_buf),
                         ImGuiInputTextFlags_CharsDecimal);
        ImGui::PopStyleColor();
        if (ImGui::IsItemDeactivatedAfterEdit() && coord_draft_buf[0] != '\0') {
          float val = static_cast<float>(expr::eval(coord_draft_buf));
          coord_draft_buf[0] = '\0';
          if (!std::isnan(val)) {
            geo.nodes.push_back({val});
            geo.domains.push_back({});
            geo.built = false;
            ApplyNodesToLengths(s);
          }
        }
      }

      ImGui::EndTable();
    }

  } else {
    char hint[64]; snprintf(hint, sizeof(hint),
        "Domain lengths  [%s]  (clear a row to remove it)", kUnits[s.geo_unit]);
    ImGui::TextDisabled("%s", hint);
    ImGui::Spacing();

    if (ImGui::BeginTable("len_tbl", 1, kTblFlags)) {
      int delete_d = -1;

      for (int d = 0; d < (int)s.geo_lengths.size(); ++d) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::SetNextItemWidth(-1.0f);
        char id[32]; snprintf(id, sizeof(id), "##ln%d", d);
        bool cleared = false;
        if (ExprInputFloat(id, &s.geo_lengths[d], &cleared)) {
          if (s.geo_lengths[d] <= 0.0f) s.geo_lengths[d] = 0.001f;
          geo.built = false;
        }
        if (cleared && (int)s.geo_lengths.size() > 1) delete_d = d;
      }

      if (delete_d >= 0) {
        s.geo_lengths.erase(s.geo_lengths.begin() + delete_d);
        if (delete_d < (int)geo.domains.size())
          geo.domains.erase(geo.domains.begin() + delete_d);
        geo.built = false;
        ApplyLengthsToNodes(s);
      }

      {
        static char len_draft_buf[32] = "";
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::SetNextItemWidth(-1.0f);
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyle().Colors[ImGuiCol_TextDisabled]);
        ImGui::InputText("##ldraft", len_draft_buf, sizeof(len_draft_buf),
                         ImGuiInputTextFlags_CharsDecimal);
        ImGui::PopStyleColor();
        if (ImGui::IsItemDeactivatedAfterEdit() && len_draft_buf[0] != '\0') {
          float val = static_cast<float>(expr::eval(len_draft_buf));
          len_draft_buf[0] = '\0';
          if (!std::isnan(val) && val > 0.0f) {
            s.geo_lengths.push_back(val);
            geo.domains.push_back({});
            geo.built = false;
          }
        }
      }

      ImGui::EndTable();
    }
  }

  ImGui::Spacing();

  bool sorted = true;
  if (s.geo_input_mode == 0) {
    for (int i = 1; i < (int)geo.nodes.size(); ++i)
      if (geo.nodes[i].x <= geo.nodes[i - 1].x) { sorted = false; break; }
  } else {
    for (float L : s.geo_lengths)
      if (L <= 0.0f) { sorted = false; break; }
  }
  bool can_build = ((int)geo.nodes.size() >= 2 || (int)s.geo_lengths.size() >= 1) && sorted;

  if (!can_build) ImGui::BeginDisabled();
  if (ImGui::Button("Build", ImVec2(-1.0f, 28))) {
    if (s.geo_input_mode == 1) ApplyLengthsToNodes(s);
    geo.built = true;
    geo.selected_node = geo.selected_domain = -1;
    s.has_results = false;
    s.status_msg.clear();
  }
  if (!can_build) ImGui::EndDisabled();

  if (!sorted) {
    ImGui::Spacing();
    ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.3f, 1.0f),
        s.geo_input_mode == 0 ? "Coordinates must be strictly increasing"
                               : "All lengths must be positive");
  }
}

// ============================================================
// Geometry View — CAD canvas
// ============================================================

static void ShowGeometryView(AppState& s) {
  static constexpr const char* kUnits[] = {"m", "cm", "mm", "\xc2\xb5m"};
  Geometry1D& geo = s.geo;

  ImDrawList* dl   = ImGui::GetWindowDrawList();
  ImVec2      orig = ImGui::GetWindowPos();
  float       cw   = ImGui::GetWindowWidth();
  float       ch   = ImGui::GetWindowHeight();

  if (s.physics == PhysicsModule::SingleParticle) {
    // Schematic of the two representative particles, sized by relative radius.
    float cy   = orig.y + ch * 0.5f;
    float rmax = std::max(s.spm.R_n, s.spm.R_p);
    if (rmax <= 0.0f) rmax = 1.0f;
    float base = std::min(cw, ch) * 0.18f;
    float rn   = base * (s.spm.R_n / rmax);
    float rp   = base * (s.spm.R_p / rmax);
    float cxn  = orig.x + cw * 0.33f;
    float cxp  = orig.x + cw * 0.67f;

    auto draw_particle = [&](float cx, float r, ImU32 col, const char* tag,
                             float radius_m) {
      for (int k = 4; k >= 1; --k) {  // concentric shells → diffusion hint
        float rr = r * k / 4.0f;
        dl->AddCircle(ImVec2(cx, cy), rr,
                      IM_COL32(((col >> 0) & 0xFF), ((col >> 8) & 0xFF),
                               ((col >> 16) & 0xFF), 60), 0, 1.0f);
      }
      dl->AddCircle(ImVec2(cx, cy), r, col, 0, 2.5f);
      dl->AddCircleFilled(ImVec2(cx, cy), 2.5f, col);
      char lbl[64];
      snprintf(lbl, sizeof(lbl), "%s   R = %.3g m", tag, (double)radius_m);
      ImVec2 ts = ImGui::CalcTextSize(lbl);
      dl->AddText(ImVec2(cx - ts.x * 0.5f, cy + r + 10.0f),
                  IM_COL32(190, 205, 230, 220), lbl);
    };

    draw_particle(cxn, rn, IM_COL32(120, 185, 255, 255), "Negative", s.spm.R_n);
    draw_particle(cxp, rp, IM_COL32(255, 170, 90, 255), "Positive", s.spm.R_p);

    const char* title = "Single Particle Model";
    ImVec2 tts = ImGui::CalcTextSize(title);
    dl->AddText(ImVec2(orig.x + (cw - tts.x) * 0.5f, orig.y + 16.0f),
                IM_COL32(150, 170, 200, 200), title);
    return;
  }

  if (s.physics == PhysicsModule::SingleParticleElectrolyte) {
    // Cell sandwich (negative | separator | positive) with the two
    // representative particles drawn above their electrodes.
    const SpmeInputs& m = s.spme;
    float Ln = m.core.L_n, Ls = m.L_s, Lp = m.core.L_p;
    float Ltot = std::max(Ln + Ls + Lp, 1e-12f);

    float pad   = 60.0f;
    float usable = cw - 2.0f * pad;
    float x0px  = orig.x + pad;
    float bandTop = orig.y + ch * 0.56f;
    float bandBot = orig.y + ch * 0.78f;
    float xn = x0px + usable * (Ln / Ltot);
    float xs = x0px + usable * ((Ln + Ls) / Ltot);
    float xe = x0px + usable;

    auto band = [&](float xa, float xb, ImU32 col, const char* lbl) {
      dl->AddRectFilled(ImVec2(xa, bandTop), ImVec2(xb, bandBot), col, 3.0f);
      dl->AddRect(ImVec2(xa, bandTop), ImVec2(xb, bandBot),
                  IM_COL32(200, 210, 230, 120), 3.0f);
      ImVec2 ts = ImGui::CalcTextSize(lbl);
      if (xb - xa > ts.x + 6.0f)
        dl->AddText(ImVec2((xa + xb) * 0.5f - ts.x * 0.5f,
                           (bandTop + bandBot) * 0.5f - ts.y * 0.5f),
                    IM_COL32(20, 25, 35, 230), lbl);
    };
    band(x0px, xn, IM_COL32(120, 185, 255, 200), "Negative");
    band(xn,   xs, IM_COL32(170, 180, 200, 180), "Separator");
    band(xs,   xe, IM_COL32(255, 170, 90, 200),  "Positive");

    // Particles above each electrode.
    float rmax = std::max(m.core.R_n, m.core.R_p);
    if (rmax <= 0.0f) rmax = 1.0f;
    float base = std::min(cw, ch) * 0.11f;
    float pcy  = orig.y + ch * 0.30f;
    auto particle = [&](float cx, float radius_m, ImU32 col) {
      float r = base * (radius_m / rmax);
      for (int k = 3; k >= 1; --k)
        dl->AddCircle(ImVec2(cx, pcy), r * k / 3.0f,
                      IM_COL32(((col >> 0) & 0xFF), ((col >> 8) & 0xFF),
                               ((col >> 16) & 0xFF), 60), 0, 1.0f);
      dl->AddCircle(ImVec2(cx, pcy), r, col, 0, 2.5f);
      dl->AddCircleFilled(ImVec2(cx, pcy), 2.0f, col);
    };
    particle((x0px + xn) * 0.5f, m.core.R_n, IM_COL32(120, 185, 255, 255));
    particle((xs + xe) * 0.5f,   m.core.R_p, IM_COL32(255, 170, 90, 255));

    dl->AddText(ImVec2(x0px, bandBot + 8.0f), IM_COL32(150, 170, 200, 200),
                "x = 0");
    const char* xl = "x = L";
    ImVec2 xls = ImGui::CalcTextSize(xl);
    dl->AddText(ImVec2(xe - xls.x, bandBot + 8.0f),
                IM_COL32(150, 170, 200, 200), xl);

    const char* title = "Single Particle Model w/ Electrolyte";
    ImVec2 tts = ImGui::CalcTextSize(title);
    dl->AddText(ImVec2(orig.x + (cw - tts.x) * 0.5f, orig.y + 16.0f),
                IM_COL32(150, 170, 200, 200), title);
    return;
  }

  if (!geo.built || (int)geo.nodes.size() < 2) {
    const char* msg = "Configure geometry and press  Build";
    ImVec2 ts = ImGui::CalcTextSize(msg);
    dl->AddText(ImVec2(orig.x + (cw - ts.x) * 0.5f, orig.y + (ch - ts.y) * 0.5f),
                IM_COL32(100, 100, 120, 200), msg);
    return;
  }

  const float pad    = 60.0f;
  const float usable = cw - 2.0f * pad;
  const float line_y = orig.y + ch * 0.50f;

  float x_min = geo.nodes.front().x;
  float x_max = geo.nodes.back().x;
  float span  = x_max - x_min;
  if (span <= 0.0f) span = 1.0f;

  auto node_px = [&](float x) -> float {
    return orig.x + pad + (x - x_min) / span * usable;
  };

  ImVec2 mouse     = ImGui::GetMousePos();
  int hovered_node   = -1;
  int hovered_domain = -1;
  bool canvas_hov  = ImGui::IsWindowHovered(ImGuiHoveredFlags_ChildWindows);

  if (canvas_hov) {
    for (int i = 0; i < (int)geo.nodes.size(); ++i) {
      float px = node_px(geo.nodes[i].x);
      float dx = mouse.x - px, dy = mouse.y - line_y;
      if (dx * dx + dy * dy < 9.0f * 9.0f) { hovered_node = i; break; }
    }
    if (hovered_node < 0) {
      for (int d = 0; d < (int)geo.domains.size(); ++d) {
        float px_l = node_px(geo.nodes[d].x);
        float px_r = node_px(geo.nodes[d + 1].x);
        if (mouse.x > px_l && mouse.x < px_r && fabsf(mouse.y - line_y) < 7.0f)
          hovered_domain = d;
      }
    }
  }

  if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
    if (hovered_node >= 0) {
      geo.selected_node = hovered_node; geo.selected_domain = -1;
    } else if (hovered_domain >= 0) {
      geo.selected_domain = hovered_domain; geo.selected_node = -1;
    } else if (canvas_hov) {
      geo.selected_node = geo.selected_domain = -1;
    }
  }

  if (hovered_node >= 0 || hovered_domain >= 0)
    ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);

  // Faint axis
  dl->AddLine(ImVec2(orig.x + pad * 0.5f, line_y),
              ImVec2(orig.x + cw - pad * 0.5f, line_y),
              IM_COL32(80, 80, 100, 60), 1.0f);

  // Domain segments
  for (int d = 0; d < (int)geo.domains.size(); ++d) {
    float px_l = node_px(geo.nodes[d].x);
    float px_r = node_px(geo.nodes[d + 1].x);
    bool  sel  = (d == geo.selected_domain);
    bool  hov  = (d == hovered_domain);
    ImU32 col  = sel ? IM_COL32(255, 190, 50, 255) :
                 hov ? IM_COL32(120, 185, 255, 255) :
                       IM_COL32(75,  145, 245, 255);
    dl->AddLine(ImVec2(px_l, line_y), ImVec2(px_r, line_y), col, 3.5f);

    char lbl[8]; snprintf(lbl, sizeof(lbl), "%d", d + 1);
    float mid_x = (px_l + px_r) * 0.5f;
    ImVec2 ts = ImGui::CalcTextSize(lbl);
    dl->AddText(ImVec2(mid_x - ts.x * 0.5f, line_y - ts.y - 10.0f),
                sel ? IM_COL32(255, 210, 100, 200) : IM_COL32(120, 160, 220, 160), lbl);
  }

  // Nodes
  for (int i = 0; i < (int)geo.nodes.size(); ++i) {
    float px       = node_px(geo.nodes[i].x);
    bool  sel      = (i == geo.selected_node);
    bool  hov      = (i == hovered_node);
    bool  terminal = (i == 0 || i == (int)geo.nodes.size() - 1);
    float r        = terminal ? 6.5f : 4.5f;

    ImU32 fill = sel ? IM_COL32(255, 200, 60, 255) :
                 hov ? IM_COL32(220, 235, 255, 255) :
                       IM_COL32(200, 220, 255, terminal ? 255 : 180);
    dl->AddCircleFilled(ImVec2(px, line_y), r, fill);
    if (sel)
      dl->AddCircle(ImVec2(px, line_y), r + 3.5f, IM_COL32(255, 200, 60, 180), 0, 1.5f);
    else if (hov)
      dl->AddCircle(ImVec2(px, line_y), r + 2.5f, IM_COL32(200, 220, 255, 120), 0, 1.0f);

    dl->AddLine(ImVec2(px, line_y + r + 1.0f), ImVec2(px, line_y + r + 5.0f),
                IM_COL32(150, 170, 200, 140), 1.0f);

    if (terminal || sel || hov) {
      char xlbl[48]; snprintf(xlbl, sizeof(xlbl), "%.4g %s",
                               (double)geo.nodes[i].x, kUnits[s.geo_unit]);
      ImVec2 xs = ImGui::CalcTextSize(xlbl);
      ImU32  tc = sel ? IM_COL32(255, 210, 100, 220) : IM_COL32(170, 185, 210, 200);
      dl->AddText(ImVec2(px - xs.x * 0.5f, line_y + r + 7.0f), tc, xlbl);
    }
  }

  // Selection overlay
  if (geo.selected_node >= 0 || geo.selected_domain >= 0) {
    const char* u = kUnits[s.geo_unit];
    char info[128] = {};
    if (geo.selected_node >= 0) {
      int  i    = geo.selected_node;
      bool term = (i == 0 || i == (int)geo.nodes.size() - 1);
      snprintf(info, sizeof(info), "Node %d  x = %.6g %s%s",
               i, (double)geo.nodes[i].x, u, term ? "" : "  (internal)");
    } else {
      int d = geo.selected_domain;
      double L = geo.nodes[d + 1].x - geo.nodes[d].x;
      snprintf(info, sizeof(info), "Domain %d  L = %.6g %s  n = %d",
               d + 1, L, u, geo.domains[d].n_cells);
    }
    ImVec2 ts = ImGui::CalcTextSize(info);
    float  ix = orig.x + 10.0f, iy = orig.y + ch - ts.y - 10.0f;
    dl->AddRectFilled(ImVec2(ix - 4, iy - 3), ImVec2(ix + ts.x + 4, iy + ts.y + 3),
                      IM_COL32(20, 20, 30, 180), 3.0f);
    dl->AddText(ImVec2(ix, iy), IM_COL32(200, 210, 230, 230), info);
  }
}

// ============================================================
// Configuration panels
// ============================================================

// Single Particle Model parameter editor.  Grouped by physical role; every
// field is editable as a number or an expression (e.g. "5e-6", "0.84*33133").
static void ShowSpmCoreInputs(SpmInputs& m, bool show_ce) {
  // Material pickers — selecting a material loads its properties (and its OCP
  // curve) into the fields below, which remain editable afterwards.
  ImGui::SeparatorText("Materials");
  ImGui::Spacing();
  static const auto kNegNames = mphys::materials::Database::ElectrodeNames(
      mphys::materials::Domain::kNegativeElectrode);
  static const auto kPosNames = mphys::materials::Database::ElectrodeNames(
      mphys::materials::Domain::kPositiveElectrode);
  ImGui::TextUnformatted("Negative electrode");
  if (MaterialCombo("##negmat", kNegNames, &m.neg_material)) ApplyNegElectrode(m);
  ImGui::Spacing();
  ImGui::TextUnformatted("Positive electrode");
  if (MaterialCombo("##posmat", kPosNames, &m.pos_material)) ApplyPosElectrode(m);
  if (show_ce) {
    static const auto kElyteNames =
        mphys::materials::Database::ElectrolyteNames();
    ImGui::Spacing();
    ImGui::TextUnformatted("Electrolyte");
    if (MaterialCombo("##elytemat", kElyteNames, &m.elyte_material)) {
      m.c_e = static_cast<float>(
          mphys::materials::Database::Electrolyte(ElectrolyteIdAt(m.elyte_material))
              .c_typical);
    }
  }

  ImGui::Spacing();
  ImGui::SeparatorText("Particle Geometry");
  ImGui::Spacing();
  LabeledFloat("Negative radius  R_n", "##Rn", &m.R_n, "m", "%.4g");
  ImGui::Spacing();
  LabeledFloat("Positive radius  R_p", "##Rp", &m.R_p, "m", "%.4g");
  ImGui::Spacing();
  LabeledFloat("Negative thickness  L_n", "##Ln", &m.L_n, "m", "%.4g");
  ImGui::Spacing();
  LabeledFloat("Positive thickness  L_p", "##Lp", &m.L_p, "m", "%.4g");
  ImGui::Spacing();
  LabeledFloat("Electrode area  A", "##area", &m.A, "m\xc2\xb2", "%.4g");
  ImGui::Spacing();
  LabeledFloat("Neg. active fraction", "##epsn", &m.eps_n, nullptr, "%.4g");
  ImGui::Spacing();
  LabeledFloat("Pos. active fraction", "##epsp", &m.eps_p, nullptr, "%.4g");

  ImGui::Spacing();
  ImGui::SeparatorText("Solid Diffusivity");
  ImGui::Spacing();
  LabeledFloat("Negative  D_n", "##Dn", &m.D_n, "m\xc2\xb2/s", "%.3e");
  ImGui::Spacing();
  LabeledFloat("Positive  D_p", "##Dp", &m.D_p, "m\xc2\xb2/s", "%.3e");

  ImGui::Spacing();
  ImGui::SeparatorText("Concentrations");
  ImGui::Spacing();
  LabeledFloat("Neg. max  c_n,max", "##cnmax", &m.cn_max, "mol/m\xc2\xb3", "%.5g");
  ImGui::Spacing();
  LabeledFloat("Pos. max  c_p,max", "##cpmax", &m.cp_max, "mol/m\xc2\xb3", "%.5g");
  ImGui::Spacing();
  LabeledFloat("Initial neg. stoich.  x0", "##x0", &m.x0, nullptr, "%.4g");
  if (m.x0 < 0.0f) m.x0 = 0.0f; if (m.x0 > 1.0f) m.x0 = 1.0f;
  ImGui::Spacing();
  LabeledFloat("Initial pos. stoich.  y0", "##y0", &m.y0, nullptr, "%.4g");
  if (m.y0 < 0.0f) m.y0 = 0.0f; if (m.y0 > 1.0f) m.y0 = 1.0f;
  if (show_ce) {
    ImGui::Spacing();
    LabeledFloat("Electrolyte conc.  c_e", "##ce", &m.c_e, "mol/m\xc2\xb3", "%.5g");
  }

  ImGui::Spacing();
  ImGui::SeparatorText("Kinetics");
  ImGui::Spacing();
  LabeledFloat("Neg. rate constant  k_n", "##kn", &m.k_n, nullptr, "%.3e");
  ImGui::Spacing();
  LabeledFloat("Pos. rate constant  k_p", "##kp", &m.k_p, nullptr, "%.3e");

  ImGui::Spacing();
  ImGui::SeparatorText("Operating Conditions");
  ImGui::Spacing();
  LabeledFloat("Applied current  I", "##cur", &m.I, "A", "%.4g");
  ImGui::SameLine();
  ImGui::TextDisabled("(discharge +)");
  ImGui::Spacing();
  LabeledFloat("Temperature  T", "##temp", &m.T, "K", "%.5g");
}

static void ShowSpmPhysics(AppState& s) {
  ImGui::TextDisabled("Each electrode is one spherical particle (PyBaMM SPM).");
  ImGui::Spacing();
  ShowSpmCoreInputs(s.spm, /*show_ce=*/true);
}

// SPMe parameter editor: the SPM particle core plus the electrolyte transport
// and separator parameters that distinguish the SPMe (PyBaMM SPMe).
static void ShowSpmePhysics(AppState& s) {
  SpmeInputs& m = s.spme;

  ImGui::TextDisabled("Two particles + electrolyte transport (PyBaMM SPMe).");
  ImGui::Spacing();
  ShowSpmCoreInputs(m.core, /*show_ce=*/false);  // electrolyte is now a field

  ImGui::Spacing();
  ImGui::SeparatorText("Separator");
  ImGui::Spacing();
  LabeledFloat("Separator thickness  L_s", "##Ls", &m.L_s, "m", "%.4g");

  ImGui::Spacing();
  ImGui::SeparatorText("Electrolyte Porosity");
  ImGui::Spacing();
  LabeledFloat("Neg. electrolyte frac.", "##een", &m.eps_e_n, nullptr, "%.4g");
  ImGui::Spacing();
  LabeledFloat("Sep. electrolyte frac.", "##ees", &m.eps_e_s, nullptr, "%.4g");
  ImGui::Spacing();
  LabeledFloat("Pos. electrolyte frac.", "##eep", &m.eps_e_p, nullptr, "%.4g");

  ImGui::Spacing();
  ImGui::SeparatorText("Electrolyte Transport");
  ImGui::Spacing();
  // Picking an electrolyte samples its concentration-dependent transport
  // functions at the nominal salt concentration to fill the constant-property
  // SPMe fields below.
  static const auto kElyteNames = mphys::materials::Database::ElectrolyteNames();
  ImGui::TextUnformatted("Electrolyte");
  if (MaterialCombo("##spmeelyte", kElyteNames, &m.core.elyte_material)) {
    const auto& e = mphys::materials::Database::Electrolyte(
        ElectrolyteIdAt(m.core.elyte_material));
    const double ce = e.c_typical;
    m.ce0     = static_cast<float>(ce);
    m.core.c_e = static_cast<float>(ce);
    m.D_e     = static_cast<float>(e.diffusivity(ce, m.core.T));
    m.kappa_e = static_cast<float>(e.conductivity(ce, m.core.T));
    m.t_plus  = static_cast<float>(e.transference(ce, m.core.T));
  }
  ImGui::Spacing();
  LabeledFloat("Salt diffusivity  D_e", "##De", &m.D_e, "m\xc2\xb2/s", "%.3e");
  ImGui::Spacing();
  LabeledFloat("Ionic conductivity  \xce\xba_e", "##kappae", &m.kappa_e, "S/m", "%.4g");
  ImGui::Spacing();
  LabeledFloat("Transference  t+", "##tplus", &m.t_plus, nullptr, "%.4g");
  ImGui::Spacing();
  LabeledFloat("Bruggeman exponent", "##brugg", &m.brugg, nullptr, "%.4g");
  ImGui::Spacing();
  LabeledFloat("Initial electrolyte conc.", "##ce0", &m.ce0, "mol/m\xc2\xb3", "%.5g");

  ImGui::Spacing();
  ImGui::SeparatorText("Solid Conductivity");
  ImGui::Spacing();
  LabeledFloat("Negative  \xcf\x83_n", "##sign", &m.sigma_n, "S/m", "%.4g");
  ImGui::Spacing();
  LabeledFloat("Positive  \xcf\x83_p", "##sigp", &m.sigma_p, "S/m", "%.4g");
}

static void ShowPhysicsPanel(AppState& s) {
  static constexpr const char* kMods[] = {
      "Convection-Diffusion-Reaction",
      "Steady Diffusion",
      "Darcy Packed Bed",
      "Single Particle Model",
      "Single Particle Model w/ Electrolyte"
  };
  bool is_darcy     = (s.physics == PhysicsModule::DarcyPackedBed);
  bool is_spm       = (s.physics == PhysicsModule::SingleParticle);
  bool is_spme      = (s.physics == PhysicsModule::SingleParticleElectrolyte);
  bool is_transient = (s.physics == PhysicsModule::ConvectionDiffusion || is_darcy);

  ImGui::SeparatorText("Physics Module");
  ImGui::Spacing();
  int sel = static_cast<int>(s.physics);
  ImGui::SetNextItemWidth(-1.0f);
  if (ImGui::Combo("##phys", &sel, kMods, IM_ARRAYSIZE(kMods))) {
    PhysicsModule prev = s.physics;
    s.physics = static_cast<PhysicsModule>(sel);
    // First switch into a particle model: seed discharge-friendly time settings.
    const bool now_particle = s.physics == PhysicsModule::SingleParticle ||
                              s.physics == PhysicsModule::SingleParticleElectrolyte;
    const bool was_particle = prev == PhysicsModule::SingleParticle ||
                              prev == PhysicsModule::SingleParticleElectrolyte;
    if (now_particle && !was_particle) {
      s.t_end = 3600.0f; s.dt_max = 20.0f; s.dt_snapshot = 20.0f;
    }
    is_spm  = (s.physics == PhysicsModule::SingleParticle);
    is_spme = (s.physics == PhysicsModule::SingleParticleElectrolyte);
  }

  ImGui::Spacing();

  if (is_spm)  { ShowSpmPhysics(s);  return; }
  if (is_spme) { ShowSpmePhysics(s); return; }

  if (!s.geo.built || s.geo.domains.empty()) {
    ImGui::TextDisabled("Build the geometry to configure physics parameters.");
  } else {
    for (int d = 0; d < (int)s.geo.domains.size(); ++d) {
      GeoDomain& dom = s.geo.domains[d];
      char sec[32]; snprintf(sec, sizeof(sec), "Domain %d", d + 1);
      ImGui::SeparatorText(sec);
      ImGui::Spacing();

      if (is_darcy) {
        char id_k[16], id_m[16];
        snprintf(id_k, sizeof(id_k), "##kp%d", d);
        snprintf(id_m, sizeof(id_m), "##mu%d", d);
        LabeledFloat("Permeability",      id_k, &dom.kappa, "m\xc2\xb2");
        ImGui::Spacing();
        LabeledFloat("Dynamic viscosity", id_m, &dom.mu,    "Pa\xc2\xb7s");
      } else {
        char id_D[16], id_u[16], id_k[16];
        snprintf(id_D, sizeof(id_D), "##D%d", d);
        snprintf(id_u, sizeof(id_u), "##u%d", d);
        snprintf(id_k, sizeof(id_k), "##k%d", d);
        LabeledFloat("Diffusivity", id_D, &dom.D, "m\xc2\xb2/s");
        if (is_transient) {
          ImGui::Spacing();
          LabeledFloat("Velocity",      id_u, &dom.u, "m/s");
          ImGui::Spacing();
          LabeledFloat("Reaction rate", id_k, &dom.k, "1/s");
        }
      }
      ImGui::Spacing();
    }
  }

  ImGui::Spacing();
  ImGui::SeparatorText("Boundary Conditions");
  ImGui::Spacing();

  if (!s.geo.built || s.geo.nodes.size() < 2) {
    ImGui::TextDisabled("Build the geometry to configure boundary conditions.");
    return;
  }

  if (is_darcy) {
    // Packed bed: each boundary is either a Pressure or Velocity specification
    auto ShowPbBc = [](const char* label, float x_pos, PackedBedBc& bc) {
      char buf[64]; snprintf(buf, sizeof(buf), "%s  (x = %.4g)", label, (double)x_pos);
      ImGui::TextUnformatted(buf);
      ImGui::Spacing();
      char r0[64], r1[64], vid[64];
      snprintf(r0,  sizeof(r0),  "Pressure##t%s",  label);
      snprintf(r1,  sizeof(r1),  "Velocity##t%s",  label);
      snprintf(vid, sizeof(vid), "##pbv%s",         label);
      ImGui::RadioButton(r0, &bc.type, 0); ImGui::SameLine();
      ImGui::RadioButton(r1, &bc.type, 1);
      ImGui::SetNextItemWidth(-1.0f - kUnitColWidth);
      float avail = ImGui::GetContentRegionAvail().x;
      const char* unit = (bc.type == 0) ? "Pa" : "m/s";
      ImGui::SetNextItemWidth(avail - ImGui::CalcTextSize(unit).x - 8.0f);
      ExprInputFloat(vid, &bc.value);
      ImGui::SameLine(0, 6); ImGui::TextDisabled("%s", unit);
    };

    ShowPbBc("Boundary 1", s.geo.nodes.front().x, s.pb_left_bc);
    ImGui::Spacing(); ImGui::Separator(); ImGui::Spacing();
    ShowPbBc("Boundary 2", s.geo.nodes.back().x,  s.pb_right_bc);
  } else {
    auto ShowBc = [](const char* label, float x_pos, NodeBc& bc) {
      char buf[64]; snprintf(buf, sizeof(buf), "%s  (x = %.4g)", label, (double)x_pos);
      ImGui::TextUnformatted(buf);
      ImGui::Spacing();
      char rid[64], nid[64], vid[64];
      snprintf(rid, sizeof(rid), "Dirichlet##t%s", label);
      snprintf(nid, sizeof(nid), "Neumann##t%s",   label);
      ImGui::RadioButton(rid, &bc.type, 0); ImGui::SameLine();
      ImGui::RadioButton(nid, &bc.type, 1);
      snprintf(vid, sizeof(vid), "##v%s", label);
      ImGui::SetNextItemWidth(-1.0f);
      ImGui::InputFloat(vid, &bc.value, 0, 0, "%.4g");
    };

    ShowBc("Boundary 1", s.geo.nodes.front().x, s.left_bc);
    ImGui::Spacing(); ImGui::Separator(); ImGui::Spacing();
    ShowBc("Boundary 2", s.geo.nodes.back().x,  s.right_bc);
  }

  if ((int)s.geo.nodes.size() > 2) {
    ImGui::Spacing();
    ImGui::TextDisabled("%d internal node(s) — continuity applied automatically.",
                        (int)s.geo.nodes.size() - 2);
  }
}

static void ShowMeshPanel(AppState& s) {
  ImGui::SeparatorText("Discretisation");
  ImGui::Spacing();

  if (s.physics == PhysicsModule::SingleParticle) {
    ImGui::TextDisabled("Radial points per particle (r/R in [0,1]).");
    ImGui::Spacing();
    LabeledInt("Cells per particle", "##spmnc", &s.spm.n_cells);
    if (s.spm.n_cells < 2) s.spm.n_cells = 2;
    ImGui::Spacing();
    ImGui::Text("Total unknowns: %d  (2 particles + voltage)",
                2 * s.spm.n_cells + 1);
    return;
  }

  if (s.physics == PhysicsModule::SingleParticleElectrolyte) {
    ImGui::TextDisabled("Electrolyte cells per region (x across the sandwich).");
    ImGui::Spacing();
    LabeledInt("Negative electrode", "##spmenn", &s.spme.n_n);
    if (s.spme.n_n < 2) s.spme.n_n = 2;
    ImGui::Spacing();
    LabeledInt("Separator", "##spmens", &s.spme.n_s);
    if (s.spme.n_s < 1) s.spme.n_s = 1;
    ImGui::Spacing();
    LabeledInt("Positive electrode", "##spmenp", &s.spme.n_p);
    if (s.spme.n_p < 2) s.spme.n_p = 2;
    ImGui::Spacing();
    const int ntot = s.spme.n_n + s.spme.n_s + s.spme.n_p;
    ImGui::TextDisabled("Particles share the total cell count as radial points.");
    ImGui::Text("Total unknowns: %d  (2 particles + electrolyte + voltage)",
                3 * ntot + 1);
    return;
  }

  if (!s.geo.built || s.geo.domains.empty()) {
    ImGui::TextDisabled("Build the geometry first (Geometry -> Build).");
    return;
  }

  if ((int)s.geo.domains.size() == 1) {
    LabeledInt("Number of cells", "##nc0", &s.geo.domains[0].n_cells);
    if (s.geo.domains[0].n_cells < 1) s.geo.domains[0].n_cells = 1;
  } else {
    if (ImGui::BeginTable("mesh_tbl", 3,
        ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
        ImGuiTableFlags_SizingStretchSame)) {
      ImGui::TableSetupColumn("Domain");
      ImGui::TableSetupColumn("Length");
      ImGui::TableSetupColumn("Cells");
      ImGui::TableHeadersRow();
      for (int d = 0; d < (int)s.geo.domains.size(); ++d) {
        double len = s.geo.nodes[d + 1].x - s.geo.nodes[d].x;
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0); ImGui::Text("%d", d + 1);
        ImGui::TableSetColumnIndex(1); ImGui::Text("%.4g", len);
        ImGui::TableSetColumnIndex(2);
        ImGui::SetNextItemWidth(-1.0f);
        char id[32]; snprintf(id, sizeof(id), "##mc%d", d);
        ImGui::InputInt(id, &s.geo.domains[d].n_cells, 0, 0);
        if (s.geo.domains[d].n_cells < 1) s.geo.domains[d].n_cells = 1;
      }
      ImGui::EndTable();
    }
  }

  ImGui::Spacing();
  int total = 0;
  for (const auto& d : s.geo.domains) total += d.n_cells;
  ImGui::Text("Total cells: %d", total);
}

static void ShowStudyPanel(AppState& s) {
  bool transient = (s.physics != PhysicsModule::SteadyDiffusion);

  ImGui::SeparatorText("Solver");
  ImGui::TextUnformatted(transient ? "Transient  (SUNDIALS IDA)"
                                   : "Steady-state  (SUNDIALS KINSOL)");

  if (transient) {
    ImGui::Spacing();
    ImGui::SeparatorText("Time Settings");
    LabeledFloat("End time",           "##tend",   &s.t_end,        "s");
    ImGui::Spacing();
    LabeledFloat("Snapshot interval",  "##dtsnap", &s.dt_snapshot,  "s");
    if (s.dt_snapshot < s.dt_max) {
      s.dt_snapshot = s.dt_max;
      ImGui::SameLine(0, 6);
      ImGui::TextDisabled("(clamped to max dt)");
    }
    ImGui::Spacing();
    ImGui::SeparatorText("Solver Settings");
    LabeledFloat("Initial dt",  "##dtinit", &s.dt_initial, "s", "%.2e");
    ImGui::Spacing();
    LabeledFloat("Max dt",      "##dtmax",  &s.dt_max,     "s");
    if (s.dt_max > s.dt_snapshot) s.dt_snapshot = s.dt_max;
  }

  ImGui::Spacing();
  ImGui::SeparatorText("Tolerances");
  LabeledFloat("Relative tolerance", "##rtol", &s.rel_tol, nullptr, "%.2e");
  ImGui::Spacing();
  LabeledFloat("Absolute tolerance", "##atol", &s.abs_tol, nullptr, "%.2e");
  ImGui::Spacing();
  ImGui::Spacing();

  bool can_run = !s.status_msg.starts_with("Running");
  if (!can_run) ImGui::BeginDisabled();
  if (ImGui::Button("Run", ImVec2(120, 32))) RunSimulation(s);
  if (!can_run) ImGui::EndDisabled();

  if (!s.status_msg.empty()) {
    ImGui::SameLine();
    ImGui::TextUnformatted(s.status_msg.c_str());
  }
}

// Linearly interpolate a field at time t between stored snapshots.
static std::vector<double> InterpolateField(const mphys::SimResult& res,
                                             int field_idx, double t) {
  const auto& snaps = res.snapshots;
  if (snaps.empty()) return {};
  if (t <= snaps.front().t) return snaps.front().fields[field_idx].values;
  if (t >= snaps.back().t)  return snaps.back().fields[field_idx].values;

  int lo = 0, hi = static_cast<int>(snaps.size()) - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (snaps[mid].t <= t) lo = mid; else hi = mid;
  }
  double t0 = snaps[lo].t, t1 = snaps[hi].t;
  double alpha = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0;
  const auto& f0 = snaps[lo].fields[field_idx].values;
  const auto& f1 = snaps[hi].fields[field_idx].values;
  std::vector<double> out(f0.size());
  for (int i = 0; i < (int)f0.size(); ++i)
    out[i] = f0[i] * (1.0 - alpha) + f1[i] * alpha;
  return out;
}

static void ShowResultsPanel(AppState& s) {
  if (!s.has_results || s.result.snapshots.empty()) {
    ImGui::TextDisabled("No results yet — run a simulation from the Study node.");
    return;
  }

  const bool is_spm  = (s.physics == PhysicsModule::SingleParticle);
  const bool is_spme = (s.physics == PhysicsModule::SingleParticleElectrolyte);
  const bool is_particle = is_spm || is_spme;

  // SPM/SPMe: terminal voltage vs time, with a marker at the selected instant.
  if (is_particle && !s.spm_voltage.empty()) {
    ImGui::SeparatorText("Terminal Voltage");
    ImGui::Spacing();
    int nv = static_cast<int>(s.spm_voltage.size());
    double v_now = s.spm_voltage.back();
    {
      // Nearest voltage sample to the current plot time.
      int best = nv - 1; double bd = 1e300;
      for (int i = 0; i < nv; ++i) {
        double d = std::abs(s.spm_time[i] - (double)s.plot_time);
        if (d < bd) { bd = d; best = i; }
      }
      v_now = s.spm_voltage[best];
    }
    ImGui::Text("V(t = %.4g s) = %.4f V", (double)s.plot_time, v_now);
    ImGui::Spacing();

    ImVec2 vsz = ImGui::GetContentRegionAvail();
    vsz.y *= 0.42f;
    if (ImPlot::BeginPlot("##spm_voltage", vsz)) {
      ImPlot::SetupAxes("t  [s]", "V  [V]",
                        ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
      ImPlot::PlotLine("V", s.spm_time.data(), s.spm_voltage.data(), nv);
      double tx = (double)s.plot_time, vy = v_now;
      ImPlot::PlotScatter("##now", &tx, &vy, 1);
      ImPlot::EndPlot();
    }
    ImGui::Spacing();
    ImGui::SeparatorText(is_spme ? "Concentration Profiles"
                                 : "Particle Concentration");
    ImGui::Spacing();
  }

  const auto& snaps   = s.result.snapshots;
  int         n_snaps = static_cast<int>(snaps.size());
  bool        transient = n_snaps > 1;

  int n_fields = static_cast<int>(snaps[0].fields.size());
  static int field_idx = 0;
  if (field_idx >= n_fields) field_idx = 0;

  if (!transient) {
    ImGui::TextUnformatted("Steady-state result");
    ImGui::Spacing();
  } else {
    float t_min = static_cast<float>(snaps.front().t);
    float t_max = static_cast<float>(snaps.back().t);

    // Clamp stored plot_time into valid range after a new run
    s.plot_time = std::clamp(s.plot_time, t_min, t_max);

    float avail = ImGui::GetContentRegionAvail().x;

    // Continuous slider across the full time range
    ImGui::SetNextItemWidth(avail - kUnitColWidth);
    ImGui::SliderFloat("##tslider", &s.plot_time, t_min, t_max, "");
    ImGui::SameLine(0, 6);
    ImGui::TextDisabled("s");
    ImGui::Spacing();

    // Exact time input via expression
    float w = avail - ImGui::CalcTextSize("t =").x - 6.0f - kUnitColWidth;
    ImGui::TextUnformatted("t ="); ImGui::SameLine(0, 6);
    ImGui::SetNextItemWidth(w);
    ExprInputFloat("##texact", &s.plot_time);
    s.plot_time = std::clamp(s.plot_time, t_min, t_max);
    ImGui::SameLine(0, 6); ImGui::TextDisabled("s");

    ImGui::Spacing();
    ImGui::TextDisabled("%d snapshots  [%.4g s … %.4g s]",
                        n_snaps, (double)t_min, (double)t_max);
    ImGui::Spacing();
  }

  if (n_fields > 1) {
    std::vector<const char*> names;
    for (const auto& f : snaps[0].fields) names.push_back(f.name.c_str());
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::Combo("##field", &field_idx, names.data(), n_fields);
    ImGui::Spacing();
  }

  // Compute values to plot — interpolated for transient, direct for steady
  std::vector<double> plot_vals;
  if (transient) {
    plot_vals = InterpolateField(s.result, field_idx, static_cast<double>(s.plot_time));
  } else {
    plot_vals = snaps[0].fields[field_idx].values;
  }

  const std::string& fname = snaps[0].fields[field_idx].name;

  // For SPMe the electrolyte field lives on the Cartesian sandwich (x in m),
  // while the particle fields live on the normalised radius r/R.  Pick the
  // matching x-axis array and label per selected field.
  const bool is_electrolyte = is_spme && fname == "c_e";
  const std::vector<double>& x_axis =
      (is_electrolyte && !s.spme_ce_x.empty()) ? s.spme_ce_x : s.cell_centres;
  const char* x_label = is_electrolyte ? "x  [m]"
                      : is_particle    ? "r / R  [-]"
                                       : "x  [m]";
  int n_cells = std::min(static_cast<int>(x_axis.size()),
                         static_cast<int>(plot_vals.size()));

  ImVec2 plot_size = ImGui::GetContentRegionAvail();
  plot_size.y -= 4.0f;

  if (ImPlot::BeginPlot("##results", plot_size)) {
    ImPlot::SetupAxes(x_label, fname.c_str(),
                      ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
    ImPlot::PlotLine(fname.c_str(),
        x_axis.data(), plot_vals.data(), n_cells);
    ImPlot::EndPlot();
  }
}

static void ShowConfigPanel(AppState& s) {
  switch (s.nav) {
    case NavNode::Geometry: ShowGeometryPanel(s); break;
    case NavNode::Physics:  ShowPhysicsPanel(s);  break;
    case NavNode::Mesh:     ShowMeshPanel(s);      break;
    case NavNode::Study:    ShowStudyPanel(s);     break;
    case NavNode::Results:  ShowResultsPanel(s);   break;
  }
}

// ============================================================
// Dock layout
// ============================================================

static void SetupDockLayout(ImGuiID dockspace_id) {
  ImGui::DockBuilderRemoveNode(dockspace_id);
  ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
  ImGui::DockBuilderSetNodeSize(dockspace_id, ImGui::GetMainViewport()->WorkSize);

  ImGuiID left, right;
  ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.18f, &left, &right);

  ImGuiID config, canvas;
  ImGui::DockBuilderSplitNode(right, ImGuiDir_Left, 0.42f, &config, &canvas);

  ImGui::DockBuilderDockWindow("Model Builder", left);
  ImGui::DockBuilderDockWindow("Configuration", config);
  ImGui::DockBuilderDockWindow("Geometry View", canvas);
  ImGui::DockBuilderFinish(dockspace_id);
}

// ============================================================
// Main render function
// ============================================================

static void RenderFrame(AppState& s) {
  static bool first_frame = true;

  ImGuiViewport* vp = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(vp->WorkPos);
  ImGui::SetNextWindowSize(vp->WorkSize);
  ImGui::SetNextWindowViewport(vp->ID);

  ImGuiWindowFlags host_flags =
      ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse |
      ImGuiWindowFlags_NoResize   | ImGuiWindowFlags_NoMove     |
      ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus |
      ImGuiWindowFlags_MenuBar;

  ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
  ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
  ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
  ImGui::Begin("##host", nullptr, host_flags);
  ImGui::PopStyleVar(3);

  if (ImGui::BeginMenuBar()) {
    if (ImGui::BeginMenu("File")) {
#ifndef __EMSCRIPTEN__
      // Native file dialogs are unavailable in the browser sandbox.
      if (ImGui::MenuItem("Open...", "Ctrl+O")) {
        static const char* kFilter[] = {"*.json"};
        const char* path = tinyfd_openFileDialog(
            "Open mphys state", "", 1, kFilter, "mphys JSON (*.json)", 0);
        if (path) {
          try {
            LoadState(s, path);
            s.status_msg = std::string("Loaded: ") + path;
          } catch (const std::exception& e) {
            s.status_msg = std::string("Error: ") + e.what();
          }
        }
      }
      if (ImGui::MenuItem("Save...", "Ctrl+S")) {
        static const char* kFilter[] = {"*.json"};
        const char* path = tinyfd_saveFileDialog(
            "Save mphys state", "untitled.json", 1, kFilter, "mphys JSON (*.json)");
        if (path) {
          try {
            SaveState(s, path);
            s.status_msg = std::string("Saved: ") + path;
          } catch (const std::exception& e) {
            s.status_msg = std::string("Error: ") + e.what();
          }
        }
      }
      ImGui::Separator();
#endif
      if (ImGui::BeginMenu("Examples")) {
        // Enumerate *.json files in MPHYS_ASSETS_DIR/examples/
        static const std::string kExDir = std::string(MPHYS_ASSETS_DIR) + "/examples/";
        // Build list once — tinyfd doesn't enumerate; use a hardcoded scan via filesystem
        // For simplicity list well-known examples embedded at compile time.
        // Add new entries here when new example JSON files are added to assets/examples/.
        static const char* kExamples[] = {
            "darcy_packed_bed.json",
            "single_particle_model.json",
            "single_particle_model_electrolyte.json",
        };
        for (const char* ex : kExamples) {
          if (ImGui::MenuItem(ex)) {
            try {
              LoadState(s, kExDir + ex);
              s.status_msg = std::string("Loaded example: ") + ex;
            } catch (const std::exception& e) {
              s.status_msg = std::string("Error: ") + e.what();
            }
          }
        }
        ImGui::EndMenu();
      }
#ifndef __EMSCRIPTEN__
      ImGui::Separator();
      if (ImGui::MenuItem("Quit", "Alt+F4"))
        glfwSetWindowShouldClose(glfwGetCurrentContext(), GLFW_TRUE);
#endif
      ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Study")) {
      if (ImGui::MenuItem("Run", "F5")) RunSimulation(s);
      ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("View")) {
      if (ImGui::MenuItem("Dark Theme",  nullptr, s.dark_theme)) {
        s.dark_theme = true;  Themes::SetDarkTheme();
      }
      if (ImGui::MenuItem("Light Theme", nullptr, !s.dark_theme)) {
        s.dark_theme = false; Themes::SetLightTheme();
      }
      ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Help")) {
      if (ImGui::BeginMenu("About")) {
        ImGui::TextUnformatted("mphys GUI — COMSOL-inspired 1D physics front-end");
        ImGui::TextUnformatted("Powered by Dear ImGui + ImPlot + SUNDIALS");
        ImGui::EndMenu();
      }
      ImGui::EndMenu();
    }

    if (!s.status_msg.empty()) {
      ImGui::SetCursorPosX(ImGui::GetCursorPosX() + 8.0f);
      ImGui::TextDisabled("%s", s.status_msg.c_str());
    }
    ImGui::EndMenuBar();
  }

  ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
  ImGui::DockSpace(dockspace_id, ImVec2(0, 0), ImGuiDockNodeFlags_None);

  if (first_frame) {
    first_frame = false;
    SetupDockLayout(dockspace_id);
  }
  ImGui::End();

  // ---- Model Builder ----
  ImGui::Begin("Model Builder", nullptr, 0);
  ImGui::TextColored(ImVec4(0.4f, 0.6f, 1.0f, 1.0f), "Model Builder");
  ImGui::Separator();
  ImGui::Spacing();

  // Taller tree-node items via FramePadding
  ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(4.0f, 7.0f));

  static const char* kNodeNames[] = {"Geometry", "Physics", "Mesh", "Study", "Results"};
  for (int i = 0; i < 5; ++i) {
    NavNode node = static_cast<NavNode>(i);
    ImGuiTreeNodeFlags flags =
        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen |
        ImGuiTreeNodeFlags_SpanFullWidth | ImGuiTreeNodeFlags_FramePadding;
    if (s.nav == node) flags |= ImGuiTreeNodeFlags_Selected;
    if (i == 4 && s.has_results)
      ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.2f, 0.75f, 0.3f, 1.0f));
    ImGui::TreeNodeEx(kNodeNames[i], flags);
    if (i == 4 && s.has_results) ImGui::PopStyleColor();
    if (ImGui::IsItemClicked()) s.nav = node;
  }

  ImGui::PopStyleVar();
  ImGui::End();

  // ---- Configuration ----
  ImGui::Begin("Configuration", nullptr, 0);
  static const char* kTitles[] = {"Geometry", "Physics", "Mesh", "Study", "Results"};
  ImGui::TextColored(ImVec4(0.4f, 0.6f, 1.0f, 1.0f), "%s", kTitles[static_cast<int>(s.nav)]);
  ImGui::Separator();
  ImGui::Spacing();
  ShowConfigPanel(s);
  ImGui::End();

  // ---- Geometry View ----
  ImGui::Begin("Geometry View", nullptr,
      ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
  ShowGeometryView(s);
  ImGui::End();
}

// ============================================================
// Render loop
// ============================================================

// Bundles everything a single frame needs, so the body can be shared between the
// desktop while-loop and the Emscripten requestAnimationFrame callback.
struct LoopContext {
  GLFWwindow* window;
  AppState* state;
};

static void MainLoopStep(void* arg) {
  LoopContext* ctx = static_cast<LoopContext*>(arg);
  GLFWwindow* window = ctx->window;

  glfwPollEvents();
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  RenderFrame(*ctx->state);

  ImGui::Render();
  int fb_w, fb_h;
  glfwGetFramebufferSize(window, &fb_w, &fb_h);
  glViewport(0, 0, fb_w, fb_h);

  const ImVec4 bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
  glClearColor(bg.x, bg.y, bg.z, bg.w);
  glClear(GL_COLOR_BUFFER_BIT);
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  glfwSwapBuffers(window);
}

// ============================================================
// Entry point
// ============================================================

int main() {
  glfwSetErrorCallback([](int error, const char* desc) {
    fprintf(stderr, "GLFW error %d: %s\n", error, desc);
  });
  if (!glfwInit()) return 1;

#if defined(__EMSCRIPTEN__)
  // Emscripten serves a WebGL2 context (GLES 3.0) via its GLFW port; the desktop
  // core-profile hints below do not apply.
  const char* glsl_version = "#version 300 es";
#elif defined(__APPLE__)
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  const char* glsl_version = "#version 150";
#else
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  const char* glsl_version = "#version 130";
#endif

#ifndef __EMSCRIPTEN__
  glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);
#endif
  GLFWwindow* window = glfwCreateWindow(1400, 900, "mphys", nullptr, nullptr);
  if (!window) { glfwTerminate(); return 1; }

  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImPlot::CreateContext();

  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

  float xscale = 1.0f;
  glfwGetWindowContentScale(window, &xscale, nullptr);
  if (xscale <= 0.0f) xscale = 1.0f;
  const float font_size = std::floor(15.0f * xscale);
  io.Fonts->AddFontFromFileTTF(
      MPHYS_ASSETS_DIR "/fonts/Roboto-VariableFont_wdth,wght.ttf", font_size);
  io.FontGlobalScale = 1.0f / xscale;

  Themes::SetDarkTheme();

  ImGui_ImplGlfw_InitForOpenGL(window, true);
#ifdef __EMSCRIPTEN__
  // Keep the GLFW window / framebuffer in sync with the HTML canvas size.
  ImGui_ImplGlfw_InstallEmscriptenCallbacks(window, "#canvas");
#endif
  ImGui_ImplOpenGL3_Init(glsl_version);

#if defined(__EMSCRIPTEN__)
  // The browser owns the event loop, so the context must outlive main(): with
  // simulate_infinite_loop the stack is unwound while the callback keeps firing.
  // Heap-allocate and intentionally never free (lives for the page's lifetime).
  auto* ctx = new LoopContext{window, new AppState()};
  emscripten_set_main_loop_arg(MainLoopStep, ctx, 0, /*simulate_infinite_loop=*/true);
#else
  AppState state;
  LoopContext ctx{window, &state};
  while (!glfwWindowShouldClose(window)) {
    MainLoopStep(&ctx);
  }

  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImPlot::DestroyContext();
  ImGui::DestroyContext();
  glfwDestroyWindow(window);
  glfwTerminate();
#endif
  return 0;
}
