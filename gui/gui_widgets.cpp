#include "gui_widgets.hpp"

#include <cmath>
#include <cstdio>
#include <unordered_map>

#include "imgui.h"

#include "expr.hpp"

namespace {

// ============================================================
// Per-widget expression input state
// ============================================================
struct ExprState { char buf[128] = ""; bool editing = false; };
std::unordered_map<ImGuiID, ExprState> g_expr_states;

}  // namespace

bool ExprInputFloat(const char* id_str, float* v, bool* cleared,
                    const char* fmt) {
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

bool ExprInputInt(const char* id_str, int* v, bool* cleared) {
  float fv = static_cast<float>(*v);
  bool  cl  = false;
  bool  ch  = ExprInputFloat(id_str, &fv, &cl);
  if (cl) { if (cleared) *cleared = true; return false; }
  if (ch) { *v = static_cast<int>(std::lround(static_cast<double>(fv))); return true; }
  return false;
}

bool LabeledFloat(const char* label, const char* id, float* v,
                  const char* unit, const char* fmt) {
  ImGui::TextUnformatted(label);
  float avail = ImGui::GetContentRegionAvail().x;
  float w = unit ? avail - kUnitColWidth : -1.0f;
  ImGui::SetNextItemWidth(w);
  bool ch = ExprInputFloat(id, v, nullptr, fmt);
  if (unit) { ImGui::SameLine(0, 6); ImGui::TextDisabled("%s", unit); }
  return ch;
}

bool LabeledInt(const char* label, const char* id, int* v, const char* unit) {
  ImGui::TextUnformatted(label);
  float avail = ImGui::GetContentRegionAvail().x;
  float w = unit ? avail - kUnitColWidth : -1.0f;
  ImGui::SetNextItemWidth(w);
  bool ch = ExprInputInt(id, v, nullptr);
  if (unit) { ImGui::SameLine(0, 6); ImGui::TextDisabled("%s", unit); }
  return ch;
}

bool MaterialCombo(const char* combo_id, const std::vector<std::string>& names,
                   int* idx) {
  if (*idx < 0 || *idx >= static_cast<int>(names.size())) *idx = 0;
  std::vector<const char*> items;
  items.reserve(names.size());
  for (const auto& n : names) items.push_back(n.c_str());
  ImGui::SetNextItemWidth(-1.0f);
  return ImGui::Combo(combo_id, idx, items.data(),
                      static_cast<int>(items.size()));
}
