#pragma once

// Application state and plain-data structures for the mphys GUI.
//
// This header holds only data (no ImGui dependency): the geometry/physics
// inputs the user edits, the selected model id, and the results buffers the
// panels read back. Material-database and SPM-parameter helpers that operate on
// this data (but do not touch the UI) live here too.

#include <map>
#include <string>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/materials/database.hpp"
#include "mphys/mesh.hpp"
#include "mphys/models/model_info.hpp"   // mphys::BcChoice
#include "mphys/models/spm.hpp"
#include "mphys/models/spme.hpp"
#include "mphys/sim_result.hpp"
#include "mphys/topology.hpp"

// Stable model ids from the registry, used for GUI dispatch.
inline constexpr const char* kConvDiffId   = "transport.conv_diff_reaction";
inline constexpr const char* kSteadyDiffId = "transport.steady_diffusion";
inline constexpr const char* kDarcyId      = "fluid_flow.darcy_packed_bed";
inline constexpr const char* kSpmId        = "liion.spm";
inline constexpr const char* kSpmeId       = "liion.spme";

// ============================================================
// 1D composite geometry (legacy StateVector model path)
// ============================================================

struct GeoNode { float x = 0.0f; };

// Per-domain physics parameters, keyed by the active model's ParamSpec::key.
// Values are filled from the model schema's defaults by EnsureModelConfig().
struct GeoDomain {
  int n_cells = 20;
  std::map<std::string, float> params;
};

struct Geometry1D {
  std::vector<GeoNode>   nodes;
  std::vector<GeoDomain> domains;
  bool built          = false;
  int  selected_node   = -1;
  int  selected_domain = -1;
};

// ============================================================
// 2D / 3D structured-box geometry (face-based MeshModel path)
// ============================================================

// Single-species convection-diffusion-reaction on a structured box. The same
// MeshModel runs on the matrix-free mesh solvers in 1D/2D/3D; the GUI uses this
// path for dim >= 2. Patch order matches MakeStructuredMesh2D/3D:
//   2D = {left, right, bottom, top}
//   3D = {left, right, bottom, top, back, front}
struct BoxGeometry {
  int dim = 2;                              // 2 or 3 (dim 1 uses Geometry1D)
  float x0 = 0.0f, x1 = 1.0f;
  float y0 = 0.0f, y1 = 0.4f;
  float z0 = 0.0f, z1 = 0.4f;
  int   nx = 60, ny = 24, nz = 16;
  // Convection-diffusion-reaction coefficients.
  float vx = 1.0f, vy = 0.0f, vz = 0.0f;    // advection velocity  [m/s]
  float D  = 0.01f;                          // diffusivity         [m^2/s]
  float k  = 1.0f;                           // reaction rate       [1/s]
  // Per-patch boundary conditions (bc_type: 0 = Dirichlet, 1 = Neumann).
  int   bc_type[6]  = {0, 1, 1, 1, 1, 1};
  float bc_value[6] = {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
  bool  steady = false;                      // steady-state vs transient
};

// One recorded field over all cells of the structured mesh at time t.
struct MeshFieldSnapshot {
  double t = 0.0;
  std::vector<double> values;   // length mesh.NCells()
};

enum class NavNode { Geometry, Physics, Mesh, Study, Results };

// ============================================================
// Single Particle Model inputs
// ============================================================

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
  // Selected database materials (indices into the domain-filtered catalogues).
  int   neg_material   = 0;                 // negative electrode
  int   pos_material   = 0;                 // positive electrode
  int   elyte_material = 0;                 // electrolyte
};

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

// ============================================================
// Database material <-> input helpers (no UI dependency)
// ============================================================

// The combos store an index into the (deterministic) domain-filtered id lists;
// these map an index back to the strongly-typed identifier.
inline mphys::materials::ElectrodeId NegElectrodeId(int idx) {
  const auto ids =
      mphys::materials::Database::ElectrodeIds(mphys::materials::Domain::kNegativeElectrode);
  if (idx < 0 || idx >= static_cast<int>(ids.size())) idx = 0;
  return ids[idx];
}
inline mphys::materials::ElectrodeId PosElectrodeId(int idx) {
  const auto ids =
      mphys::materials::Database::ElectrodeIds(mphys::materials::Domain::kPositiveElectrode);
  if (idx < 0 || idx >= static_cast<int>(ids.size())) idx = 0;
  return ids[idx];
}
inline mphys::materials::ElectrolyteId ElectrolyteIdAt(int idx) {
  const auto ids = mphys::materials::Database::ElectrolyteIds();
  if (idx < 0 || idx >= static_cast<int>(ids.size())) idx = 0;
  return ids[idx];
}

// Load a material's scalar properties into the editable input fields.
inline void ApplyNegElectrode(SpmInputs& in) {
  const auto& m = mphys::materials::Database::Electrode(NegElectrodeId(in.neg_material));
  in.cn_max = static_cast<float>(m.c_max);
  in.D_n    = static_cast<float>(m.diffusivity);
  in.k_n    = static_cast<float>(m.reaction_rate);
}
inline void ApplyPosElectrode(SpmInputs& in) {
  const auto& m = mphys::materials::Database::Electrode(PosElectrodeId(in.pos_material));
  in.cp_max = static_cast<float>(m.c_max);
  in.D_p    = static_cast<float>(m.diffusivity);
  in.k_p    = static_cast<float>(m.reaction_rate);
}

inline mphys::models::SpmParameters ToSpmParameters(const SpmInputs& in) {
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
  // are curves, not editable scalars, so the database is their only source.
  p.Un = mphys::materials::Database::Electrode(NegElectrodeId(in.neg_material)).ocp;
  p.Up = mphys::materials::Database::Electrode(PosElectrodeId(in.pos_material)).ocp;
  return p;
}

inline mphys::models::SpmeParameters ToSpmeParameters(const SpmeInputs& in) {
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

// ============================================================
// Top-level application state
// ============================================================

struct AppState {
  NavNode     nav = NavNode::Geometry;
  std::string model_id = kConvDiffId;

  mphys::CoordSystem coord_system = mphys::CoordSystem::kCartesian;

  Geometry1D geo = {
    {{0.0f}, {1.0f}},
    {{20, {}}},
    false, -1, -1
  };

  int geo_input_mode = 0;                        // 0=coordinates, 1=lengths
  std::vector<float> geo_lengths = {1.0f};
  int geo_unit = 0;                              // 0=m  1=cm  2=mm  3=µm

  int  dim = 1;                                  // 1, 2 or 3 (geometry dimension)
  BoxGeometry box;                               // 2D/3D structured-box config

  // Boundary-condition choices keyed by the active model's BcSlot::key.
  std::map<std::string, mphys::BcChoice> bcs;

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

  mphys::SimResult     result;          // 1D results (Field-based snapshots)
  std::vector<double>  cell_centres;

  // 2D/3D results (face-based mesh path).
  mphys::Mesh                     mesh_result;   // structured mesh of last run
  std::vector<MeshFieldSnapshot>  mesh_snaps;    // field "c" over time
  int                             result_dim = 1;  // dim of the data in results

  float                plot_time    = 0.0f;
  bool                 has_results  = false;
  std::string          status_msg;

  bool dark_theme = true;
};
