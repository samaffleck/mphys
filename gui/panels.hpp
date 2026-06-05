#pragma once

// The COMSOL-style configuration panels and the geometry/results canvas.

#include "app_state.hpp"

// Left/centre configuration forms (panels_config.cpp).
void ShowGeometryPanel(AppState& s);
void ShowPhysicsPanel(AppState& s);
void ShowMeshPanel(AppState& s);
void ShowStudyPanel(AppState& s);

// Right-hand canvas: geometry preview and result visualisation (panels_view.cpp).
void ShowGeometryView(AppState& s);
void ShowResultsPanel(AppState& s);
