#pragma once

// Geometry/parameter preparation and the simulation runners that drive the
// SUNDIALS solvers from the GUI's AppState.

#include "app_state.hpp"
#include "mphys/models/model_info.hpp"

// Ensure every domain and boundary slot has a value for the active model's
// parameters, filling any missing entries from the schema defaults. Existing
// values are preserved. Falls back to the first registered model if model_id is
// unknown (e.g. a legacy save). Returns the resolved ModelInfo (or nullptr).
const mphys::ModelInfo* EnsureModelConfig(AppState& s);

// Convert between the per-domain length list and the node-coordinate list.
void ApplyLengthsToNodes(AppState& s);
void ApplyNodesToLengths(AppState& s);

// Run the simulation for the currently selected model and geometry. Dispatches
// to the bespoke particle runners, the 1D registry path, or the 2D/3D mesh
// path as appropriate. Results are written back into AppState.
void RunSimulation(AppState& s);
