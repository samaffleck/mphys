#pragma once

// Reusable ImGui input widgets shared across the configuration panels:
// expression-aware numeric inputs and the standardised "label above, input box
// + inline unit" layout helpers.

#include <string>
#include <vector>

// Standard width reserved for the inline unit label to the right of an input.
inline constexpr float kUnitColWidth = 44.0f;

// Expression-aware float input. When idle shows the formatted number; while
// focused lets the user type a full expression (see expr::eval). On commit,
// evaluates and updates *v. If cleared to empty, *cleared is set (value kept).
bool ExprInputFloat(const char* id_str, float* v, bool* cleared = nullptr,
                    const char* fmt = "%.6g");

bool ExprInputInt(const char* id_str, int* v, bool* cleared = nullptr);

// "label above, input + inline unit" rows. Return true if the value changed.
bool LabeledFloat(const char* label, const char* id, float* v,
                  const char* unit = nullptr, const char* fmt = "%.6g");

bool LabeledInt(const char* label, const char* id, int* v,
                const char* unit = nullptr);

// Full-width dropdown over a list of display names. Returns true if the
// selection changed (and writes the new index into *idx).
bool MaterialCombo(const char* combo_id, const std::vector<std::string>& names,
                   int* idx);
