#include "panels.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "imgui.h"

#include "expr.hpp"
#include "gui_widgets.hpp"
#include "mphys/materials/database.hpp"
#include "mphys/models/model_registry.hpp"
#include "simulation.hpp"

// ============================================================
// Shared 2D/3D structured-box helpers
// ============================================================

static int BoxPatchCount(int dim) { return dim == 3 ? 6 : 4; }

static const char* BoxPatchName(int dim, int p) {
  static const char* n2[] = {"Left  (x min)", "Right (x max)",
                             "Bottom (y min)", "Top   (y max)"};
  static const char* n3[] = {"Left  (x min)", "Right (x max)",
                             "Bottom (y min)", "Top   (y max)",
                             "Back  (z min)",  "Front (z max)"};
  return (dim == 3 ? n3 : n2)[p];
}

// Dimension selector (1D / 2D / 3D). Returns true if the dimension changed.
static bool ShowDimensionSelector(AppState& s) {
  ImGui::SeparatorText("Dimension");
  ImGui::Spacing();
  int dsel = s.dim - 1;
  bool changed = false;
  ImGui::TextUnformatted("Geometry:");
  ImGui::SameLine();
  if (ImGui::RadioButton("1D", &dsel, 0)) { s.dim = 1; changed = true; }
  ImGui::SameLine();
  if (ImGui::RadioButton("2D", &dsel, 1)) { s.dim = 2; changed = true; }
  ImGui::SameLine();
  if (ImGui::RadioButton("3D", &dsel, 2)) { s.dim = 3; changed = true; }
  if (changed) {
    if (s.dim >= 2) s.box.dim = s.dim;
    s.has_results = false;
    s.status_msg.clear();
  }
  return changed;
}

// ============================================================
// Geometry panel
// ============================================================

static void ShowBoxExtents(AppState& s) {
  BoxGeometry& b = s.box;
  ImGui::Spacing();
  ImGui::SeparatorText("Domain Extents");
  ImGui::Spacing();
  LabeledFloat("x min", "##bx0", &b.x0, "m");
  ImGui::Spacing();
  LabeledFloat("x max", "##bx1", &b.x1, "m");
  ImGui::Spacing();
  LabeledFloat("y min", "##by0", &b.y0, "m");
  ImGui::Spacing();
  LabeledFloat("y max", "##by1", &b.y1, "m");
  if (b.dim == 3) {
    ImGui::Spacing();
    LabeledFloat("z min", "##bz0", &b.z0, "m");
    ImGui::Spacing();
    LabeledFloat("z max", "##bz1", &b.z1, "m");
  }
  ImGui::Spacing();
  ImGui::Spacing();
  ImGui::PushTextWrapPos(0.0f);
  ImGui::TextDisabled(
      "Single-species convection-diffusion-reaction on a structured box. "
      "Set velocity, D, k and boundary conditions in Physics; grid resolution "
      "in Mesh; then Run from Study.");
  ImGui::PopTextWrapPos();
}

void ShowGeometryPanel(AppState& s) {
  static constexpr const char* kUnits[] = {"m", "cm", "mm", "\xc2\xb5m"};
  Geometry1D& geo = s.geo;

  if (s.model_id == kSpmId) {
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

  if (s.model_id == kSpmeId) {
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

  // The convection-diffusion-reaction model runs in 1D / 2D / 3D. Offer the
  // dimension selector; 2D/3D switch to the structured-box mesh path.
  if (s.model_id == kConvDiffId) {
    ShowDimensionSelector(s);
    if (s.dim >= 2) { ShowBoxExtents(s); return; }
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
// Physics panel
// ============================================================

// Single Particle Model parameter editor. Grouped by physical role; every
// field is editable as a number or an expression (e.g. "5e-6", "0.84*33133").
static void ShowSpmCoreInputs(SpmInputs& m, bool show_ce) {
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

// Coefficients + boundary conditions for the 2D/3D box convection-diffusion-
// reaction model.
static void ShowBoxPhysics(AppState& s) {
  BoxGeometry& b = s.box;
  ImGui::SeparatorText("Transport Coefficients");
  ImGui::Spacing();
  ImGui::TextDisabled("dc/dt + div(u c) - D lap(c) + k c = 0");
  ImGui::Spacing();
  LabeledFloat("Velocity  u_x", "##vx", &b.vx, "m/s", "%.4g");
  ImGui::Spacing();
  LabeledFloat("Velocity  u_y", "##vy", &b.vy, "m/s", "%.4g");
  if (b.dim == 3) {
    ImGui::Spacing();
    LabeledFloat("Velocity  u_z", "##vz", &b.vz, "m/s", "%.4g");
  }
  ImGui::Spacing();
  LabeledFloat("Diffusivity  D", "##bD", &b.D, "m\xc2\xb2/s", "%.4g");
  ImGui::Spacing();
  LabeledFloat("Reaction rate  k", "##bk", &b.k, "1/s", "%.4g");

  ImGui::Spacing();
  ImGui::SeparatorText("Boundary Conditions");
  ImGui::Spacing();
  const int np = BoxPatchCount(b.dim);
  for (int p = 0; p < np; ++p) {
    ImGui::TextUnformatted(BoxPatchName(b.dim, p));
    ImGui::Spacing();
    char r0[32], r1[32], vid[32];
    snprintf(r0, sizeof(r0), "Dirichlet##d%d", p);
    snprintf(r1, sizeof(r1), "Neumann##n%d", p);
    snprintf(vid, sizeof(vid), "##bcv%d", p);
    ImGui::RadioButton(r0, &b.bc_type[p], 0);
    ImGui::SameLine();
    ImGui::RadioButton(r1, &b.bc_type[p], 1);
    float avail = ImGui::GetContentRegionAvail().x;
    ImGui::SetNextItemWidth(avail - kUnitColWidth);
    ExprInputFloat(vid, &b.bc_value[p]);
    ImGui::SameLine(0, 6);
    ImGui::TextDisabled(b.bc_type[p] == 0 ? "value" : "flux");
    if (p + 1 < np) { ImGui::Spacing(); ImGui::Separator(); ImGui::Spacing(); }
  }
}

void ShowPhysicsPanel(AppState& s) {
  const mphys::ModelInfo* info = EnsureModelConfig(s);
  auto& reg = mphys::BuiltinModels();

  ImGui::SeparatorText("Physics Model");
  ImGui::Spacing();

  // Package → Model selection tree (COMSOL-style).
  ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(4.0f, 4.0f));
  for (const auto& pkg : reg.Packages()) {
    ImGuiTreeNodeFlags pkg_flags =
        ImGuiTreeNodeFlags_SpanFullWidth | ImGuiTreeNodeFlags_DefaultOpen |
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (ImGui::TreeNodeEx(pkg.c_str(), pkg_flags)) {
      for (const mphys::ModelInfo* m : reg.InPackage(pkg)) {
        ImGuiTreeNodeFlags leaf_flags =
            ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen |
            ImGuiTreeNodeFlags_SpanFullWidth;
        if (m->id == s.model_id) leaf_flags |= ImGuiTreeNodeFlags_Selected;
        char label[128];
        snprintf(label, sizeof(label), "%s##%s", m->name.c_str(), m->id.c_str());
        ImGui::TreeNodeEx(label, leaf_flags);
        if (ImGui::IsItemClicked() && m->id != s.model_id) {
          std::string prev = s.model_id;
          s.model_id = m->id;
          info = EnsureModelConfig(s);
          const bool now_particle = s.model_id == kSpmId || s.model_id == kSpmeId;
          const bool was_particle = prev == kSpmId || prev == kSpmeId;
          if (now_particle && !was_particle) {
            s.t_end = 3600.0f; s.dt_max = 20.0f; s.dt_snapshot = 20.0f;
          }
        }
      }
      ImGui::TreePop();
    }
  }
  ImGui::PopStyleVar();

  if (info && !info->description.empty()) {
    ImGui::Spacing();
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("%s", info->description.c_str());
    ImGui::PopTextWrapPos();
  }
  ImGui::Spacing();

  if (!info) {
    ImGui::TextDisabled("No physics model is registered.");
    return;
  }

  // Models with bespoke geometry/parameters render their own panel.
  if (s.model_id == kSpmId)  { ShowSpmPhysics(s);  return; }
  if (s.model_id == kSpmeId) { ShowSpmePhysics(s); return; }

  // 2D/3D convection-diffusion-reaction: structured-box coefficients + BCs.
  if (s.model_id == kConvDiffId && s.dim >= 2) { ShowBoxPhysics(s); return; }

  // ── Per-domain parameters (driven by the model schema) ─────────────────────
  if (!s.geo.built || s.geo.domains.empty()) {
    ImGui::TextDisabled("Build the geometry to configure physics parameters.");
  } else {
    for (int d = 0; d < (int)s.geo.domains.size(); ++d) {
      GeoDomain& dom = s.geo.domains[d];
      char sec[32]; snprintf(sec, sizeof(sec), "Domain %d", d + 1);
      ImGui::SeparatorText(sec);
      ImGui::Spacing();
      for (const auto& p : info->schema.params) {
        if (p.scope != mphys::ParamScope::kPerDomain) continue;
        float& val = dom.params[p.key];
        char id[64]; snprintf(id, sizeof(id), "##p_%d_%s", d, p.key.c_str());
        LabeledFloat(p.label.c_str(), id, &val,
                     p.unit.empty() ? nullptr : p.unit.c_str());
        ImGui::Spacing();
      }
    }
  }

  // ── Boundary conditions (driven by the model schema) ───────────────────────
  ImGui::Spacing();
  ImGui::SeparatorText("Boundary Conditions");
  ImGui::Spacing();

  if (!s.geo.built || s.geo.nodes.size() < 2) {
    ImGui::TextDisabled("Build the geometry to configure boundary conditions.");
    return;
  }

  for (int bi = 0; bi < (int)info->schema.boundaries.size(); ++bi) {
    const mphys::BcSlot& slot = info->schema.boundaries[bi];
    mphys::BcChoice& choice = s.bcs[slot.key];

    float x_pos = (bi == 0) ? s.geo.nodes.front().x : s.geo.nodes.back().x;
    char hdr[96]; snprintf(hdr, sizeof(hdr), "%s  (x = %.4g)",
                           slot.label.c_str(), (double)x_pos);
    ImGui::TextUnformatted(hdr);
    ImGui::Spacing();

    for (int o = 0; o < (int)slot.options.size(); ++o) {
      char rid[80]; snprintf(rid, sizeof(rid), "%s##bc_%s_%d",
                             slot.options[o].label.c_str(), slot.key.c_str(), o);
      ImGui::RadioButton(rid, &choice.option, o);
      if (o + 1 < (int)slot.options.size()) ImGui::SameLine();
    }

    if (choice.option < 0 || choice.option >= (int)slot.options.size())
      choice.option = 0;
    const std::string& unit = slot.options[choice.option].unit;
    float avail = ImGui::GetContentRegionAvail().x;
    ImGui::SetNextItemWidth(unit.empty() ? -1.0f : avail - kUnitColWidth);
    char vid[48]; snprintf(vid, sizeof(vid), "##bcv_%s", slot.key.c_str());
    float fv = static_cast<float>(choice.value);
    if (ExprInputFloat(vid, &fv)) choice.value = fv;
    if (!unit.empty()) { ImGui::SameLine(0, 6); ImGui::TextDisabled("%s", unit.c_str()); }

    if (bi + 1 < (int)info->schema.boundaries.size()) {
      ImGui::Spacing(); ImGui::Separator(); ImGui::Spacing();
    }
  }

  if ((int)s.geo.nodes.size() > 2) {
    ImGui::Spacing();
    ImGui::TextDisabled("%d internal node(s) — continuity applied automatically.",
                        (int)s.geo.nodes.size() - 2);
  }
}

// ============================================================
// Mesh panel
// ============================================================

void ShowMeshPanel(AppState& s) {
  ImGui::SeparatorText("Discretisation");
  ImGui::Spacing();

  if (s.model_id == kSpmId) {
    ImGui::TextDisabled("Radial points per particle (r/R in [0,1]).");
    ImGui::Spacing();
    LabeledInt("Cells per particle", "##spmnc", &s.spm.n_cells);
    if (s.spm.n_cells < 2) s.spm.n_cells = 2;
    ImGui::Spacing();
    ImGui::Text("Total unknowns: %d  (2 particles + voltage)",
                2 * s.spm.n_cells + 1);
    return;
  }

  if (s.model_id == kSpmeId) {
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

  // 2D/3D structured-box resolution.
  if (s.model_id == kConvDiffId && s.dim >= 2) {
    BoxGeometry& b = s.box;
    ImGui::TextDisabled("Structured grid cells per axis.");
    ImGui::Spacing();
    LabeledInt("Cells in x", "##nx", &b.nx);
    if (b.nx < 1) b.nx = 1;
    ImGui::Spacing();
    LabeledInt("Cells in y", "##ny", &b.ny);
    if (b.ny < 1) b.ny = 1;
    if (b.dim == 3) {
      ImGui::Spacing();
      LabeledInt("Cells in z", "##nz", &b.nz);
      if (b.nz < 1) b.nz = 1;
    }
    ImGui::Spacing();
    const long total =
        (long)b.nx * b.ny * (b.dim == 3 ? b.nz : 1);
    ImGui::Text("Total cells: %ld", total);
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

// ============================================================
// Study panel
// ============================================================

void ShowStudyPanel(AppState& s) {
  const mphys::ModelInfo* info = mphys::BuiltinModels().Find(s.model_id);
  const bool mesh_path = (s.model_id == kConvDiffId && s.dim >= 2);
  bool transient = !info || info->solver == mphys::SolverKind::kTransient;
  if (mesh_path) transient = !s.box.steady;

  ImGui::SeparatorText("Solver");
  if (mesh_path) {
    ImGui::TextUnformatted(transient
        ? "Transient  (SUNDIALS IDA, matrix-free)"
        : "Steady-state  (SUNDIALS KINSOL, matrix-free)");
    ImGui::Spacing();
    ImGui::Checkbox("Steady-state solve", &s.box.steady);
    transient = !s.box.steady;
  } else {
    ImGui::TextUnformatted(transient ? "Transient  (SUNDIALS IDA)"
                                     : "Steady-state  (SUNDIALS KINSOL)");
  }

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
