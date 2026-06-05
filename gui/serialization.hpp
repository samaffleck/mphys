#pragma once

// JSON save/load of the GUI application state (cereal-based).

#include <string>

#include "app_state.hpp"

void SaveState(const AppState& s, const std::string& path);
void LoadState(AppState& s, const std::string& path);
