#include "serialization.hpp"

#include <fstream>
#include <stdexcept>

#include <cereal/archives/json.hpp>
#include <cereal/types/common.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include "simulation.hpp"   // EnsureModelConfig

// ============================================================
// Cereal serialisation — non-intrusive, JSON-friendly
//
// Note: only the 1D (Geometry1D) configuration is persisted. The 2D/3D box
// configuration is session-only, so existing 1D example JSON files load
// unchanged.
// ============================================================

template<class Ar> void serialize(Ar& ar, GeoNode& v) {
  ar(cereal::make_nvp("x", v.x));
}

template<class Ar> void serialize(Ar& ar, GeoDomain& v) {
  ar(cereal::make_nvp("n_cells", v.n_cells),
     cereal::make_nvp("params",  v.params));
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

// In mphys namespace so cereal finds it via ADL on mphys::BcChoice.
namespace mphys {
template<class Ar> void serialize(Ar& ar, BcChoice& v) {
  ar(cereal::make_nvp("option", v.option),
     cereal::make_nvp("value",  v.value));
}
}  // namespace mphys

// Split save/load for AppState so the coord-system enum round-trips as an int.
template<class Ar> void save(Ar& ar, const AppState& v) {
  ar(cereal::make_nvp("model_id",       v.model_id),
     cereal::make_nvp("coord_system",   static_cast<int>(v.coord_system)),
     cereal::make_nvp("geo",            v.geo),
     cereal::make_nvp("geo_input_mode", v.geo_input_mode),
     cereal::make_nvp("geo_lengths",    v.geo_lengths),
     cereal::make_nvp("geo_unit",       v.geo_unit),
     cereal::make_nvp("bcs",            v.bcs),
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
  int coord_i = 0;
  ar(cereal::make_nvp("model_id",       v.model_id),
     cereal::make_nvp("coord_system",   coord_i),
     cereal::make_nvp("geo",            v.geo),
     cereal::make_nvp("geo_input_mode", v.geo_input_mode),
     cereal::make_nvp("geo_lengths",    v.geo_lengths),
     cereal::make_nvp("geo_unit",       v.geo_unit),
     cereal::make_nvp("bcs",            v.bcs),
     cereal::make_nvp("spm",            v.spm),
     cereal::make_nvp("spme",           v.spme),
     cereal::make_nvp("t_end",          v.t_end),
     cereal::make_nvp("dt_initial",     v.dt_initial),
     cereal::make_nvp("dt_max",         v.dt_max),
     cereal::make_nvp("dt_snapshot",    v.dt_snapshot),
     cereal::make_nvp("rel_tol",        v.rel_tol),
     cereal::make_nvp("abs_tol",        v.abs_tol),
     cereal::make_nvp("dark_theme",     v.dark_theme));
  v.coord_system = static_cast<mphys::CoordSystem>(coord_i);
}

void SaveState(const AppState& s, const std::string& path) {
  std::ofstream os(path);
  if (!os) throw std::runtime_error("Cannot open for writing: " + path);
  cereal::JSONOutputArchive ar(os);
  ar(cereal::make_nvp("mphys_state", const_cast<AppState&>(s)));
}

void LoadState(AppState& s, const std::string& path) {
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
  s.mesh_snaps.clear();
  s.result_dim  = 1;
  s.plot_time   = 0.0f;
  s.status_msg.clear();
  s.geo.selected_node = s.geo.selected_domain = -1;
  EnsureModelConfig(s);  // fill any parameters/BCs the saved file omitted
}
