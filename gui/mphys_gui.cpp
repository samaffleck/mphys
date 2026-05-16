// /**
// *** Origional Author:       Sam Affleck
// *** Created:                07.01.2025
// **/

// // AdSim GUI Includes
// #include "MenuBar.h"
// #include "Plotting.h"
// #include "Solver.h"
// #include "Fluids.h"
// #include "ReactorMenu.h"
// #include "CycleMenu.h"
// #include "PorousMediaMenu.h"
// #include "StokesianUtil.h"
// #include "WindowHeaderTitles.h"
// #include "themes.h"
// #include "NavigationTree.h"
// #include "system_info_menu.h"
// #include "kpi_config_menu.h"
// #include "optimisation_menu.h"
// #include "parameterisation_menu.h"
// #include "cost_estimator_menu.h"
// #include "ParametersMenu.h"
// #include "parametric_study_menu.h"

// // ImGui Includes
// #include "hello_imgui/hello_imgui.h"

// // STL Includes
// #include <algorithm>

// // ImPlot Includes
// #include "implot.h"
// #include "implot_internal.h"

// // StokeSolver Includes
// #include "adsim/solver/coupled_solver.h"
// #include "adsim/solver/log.h"
// #include "adsim/objects/user_data.h"

// // STL Includes
// #include <string>

// // Emscripten Includes
// #ifdef EMSCRIPTEN
// #include <emscripten.h>
// #include <emscripten/bind.h>

// EM_JS(void, notify_ready, (), {
//   try { parent.postMessage({ type: 'wasmReady' }, window.location.origin); } catch (e) {}
// });
// #endif


// static adsim::SolverRunState solver_state;
// static adsim::OptimisationRunState optimisation_state;
// static adsim::ParameterisationRunState parameterisation_state;
// static adsim::ParametricStudyRunState parametric_study_state;

// static bool IsAnyJobRunning() {
//   return solver_state.IsRunning()
//       || optimisation_state.is_running
//       || parameterisation_state.is_running
//       || parametric_study_state.is_running;
// }

// static void ShowGUI(adsim::UserData& system) {
//   static bool on_first_load = true;

//   if (on_first_load) {
//     on_first_load = false;
//     Themes::SetFluentLight();
//     ImGuiIO& io = ImGui::GetIO();
//     io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
//     #ifdef EMSCRIPTEN
//     notify_ready();
//     #endif
//   }

//   static adsim::NavigationState nav_state;

//   // Make the variable registry available to all ShowExpressionInput widgets this frame
//   adsim::SetExpressionRegistry(&system.variable_registry);

//   // Update solver log each frame (detects status transitions)
//   adsim::UpdateSolverLog(solver_state);

//   adsim::MenuItem selected_menu = adsim::ShowNavigationTree(nav_state, system, solver_state);

//   // Show solver window (with all tabs) and coupled Results window
//   if (selected_menu == adsim::MenuItem::SolverData) {
//     adsim::ShowSolverControlsWindow(system, solver_state);
//     adsim::ShowResultsWindowForSolver(solver_state, system);
//   }

//   if (selected_menu == adsim::MenuItem::Parameters) {
//     adsim::ShowParametersMenu(system);
//   }

//   bool opt_open = selected_menu == adsim::MenuItem::Optimisation;
//   if (opt_open) {
//     adsim::ShowOptimisationMenu(system, nav_state, optimisation_state);
//   }

//   bool param_open = selected_menu == adsim::MenuItem::Parameterisation;
//   if (param_open) {
//     adsim::ShowParameterisationMenu(system, nav_state, parameterisation_state);
//   }

//   bool cost_est_open = selected_menu == adsim::MenuItem::CostEstimator;
//   if (cost_est_open) {
//     adsim::ShowCostEstimatorMenu(system);
//   }

//   // bool ps_open = selected_menu == adsim::MenuItem::ParametricStudy;
//   // if (ps_open) {
//   //   adsim::ShowParametricStudyMenu(system, nav_state, parametric_study_state);
//   // }

//   ImGui::BeginDisabled(IsAnyJobRunning());

//   bool sys_info_open = (selected_menu == adsim::MenuItem::SystemInfo);
//   if (sys_info_open) {
//     adsim::ShowSystemInformation(system.info);
//   }

//   bool fluid_tree_open = (selected_menu == adsim::MenuItem::FluidPackageSelector || selected_menu == adsim::MenuItem::Component);
//   if (fluid_tree_open) {
//     adsim::ShowFluidPackageSelector(system.fluid, system.bed, system.adsorbent_directory, nav_state);
//   }

//   bool library_tree_open = (selected_menu == adsim::MenuItem::AdsorbentLibrary || selected_menu == adsim::MenuItem::AdsorbentMaterial);
//   if (library_tree_open) {
//     ImGui::Begin(WindowHeaderTitle::ADSORBENT_LIBRARY);
//     adsim::ShowWindowHeader("Adsorbent Library");
//     adsim::ShowAdsorbentLibraryMenu(system.adsorbent_directory, system.fluid, nav_state);
//     ImGui::End();
//   }

//   bool adsorbent_property_open = selected_menu == adsim::MenuItem::AdsorbentMaterial
//       && system.adsorbent_directory.Contains(nav_state.selected_adsorbent_id);
//   if (adsorbent_property_open) {
//     ImGui::Begin(WindowHeaderTitle::ADSORBENT_MATERIAL);
//     adsim::ShowWindowHeader("Adsorbent");
//     adsim::Adsorbent& selected_adsorbent = system.adsorbent_directory.Get(nav_state.selected_adsorbent_id);
//     adsim::ShowAdsorbentMaterialMenu(selected_adsorbent, system.fluid);
//     ImGui::End();
//   }

//   bool reactor_tree_open = (selected_menu == adsim::MenuItem::ReactorProperties || selected_menu == adsim::MenuItem::AdsorbentLayer);
//   if (reactor_tree_open) {
//     adsim::ShowReactorPropertiesMenu(system.bed, system.adsorbent_directory, system.fluid, nav_state);
//   }

//   bool cycle_tree_open = (selected_menu == adsim::MenuItem::CycleDesign || selected_menu == adsim::MenuItem::CycleStep);
//   if (cycle_tree_open) {
//     adsim::ShowCycleDesignMenu(system, nav_state);
//   }

//   switch (selected_menu) {
//     case adsim::MenuItem::Component:
//       // Show the selected component menu
//       if (!nav_state.selected_component_name.empty()) {
//         // Check if component still exists in the fluid
//         auto it = std::find(system.fluid.components_.begin(), 
//                            system.fluid.components_.end(), 
//                            nav_state.selected_component_name);
//         if (it != system.fluid.components_.end()) {
//           ImGui::Begin(WindowHeaderTitle::COMPONENT_PROPERTIES);
//           adsim::ShowWindowHeader("Component Properties");
//           adsim::ShowComponentMenu(system.fluid, nav_state.selected_component_name);
//           ImGui::End();
//         } else {
//           // Component was removed, clear selection
//           nav_state.selected_component_name.clear();
//           nav_state.selected_item = adsim::MenuItem::FluidPackageSelector;
//         }
//       }
//       break;
//     case adsim::MenuItem::AdsorbentLayer:
//       // Show the selected layer menu
//       if (nav_state.selected_layer_index >= 0 &&
//           nav_state.selected_layer_index < static_cast<int>(system.bed.layers_.layers_.size())) {
//         ImGui::Begin(WindowHeaderTitle::ADSORBENT_PROPERTIES);
//         adsim::ShowWindowHeader("Adsorbent Properties");
//         auto& selectedLayer = system.bed.layers_.layers_[nav_state.selected_layer_index];
//         adsim::ShowPorousMediaLayerDataMenu(selectedLayer, system.adsorbent_directory, system.fluid, system.bed.GetArea());
//         ImGui::End();
//       }
//       break;
//     case adsim::MenuItem::CycleStep:
//       // Show the selected step menu
//       if (nav_state.selected_step_index >= 0 && 
//           nav_state.selected_step_index < static_cast<int>(system.cycle.step_list_.size())) {
//         auto& selectedStep = system.cycle.step_list_[nav_state.selected_step_index];
//         ImGui::Begin(WindowHeaderTitle::STEP_DATA);
//         adsim::ShowWindowHeader("Step Data");
//         adsim::ShowStepMenu(system.cycle, selectedStep, system.fluid, system.bed);
//         ImGui::End();
//       }
//       break;
//     default:
//       break;
//   }
  
//   ImGui::EndDisabled();
  
//   // Stokesian::ShowStatusBar(system, solver);
// }

// int main(int , char *[]) {
//   std::string app_name = "Skarstrom 2026.3.0";
//   adsim::UserData system;

//   ImPlot::CreateContext();

//   HelloImGui::RunnerParams params;
//   params.callbacks.ShowGui = [&system]() { ShowGUI(system); };
//   params.callbacks.LoadAdditionalFonts = adsim::InitializeImGuiFonts;
//   params.callbacks.ShowMenus = [&system]() { adsim::ShowMainMenuBar(system, IsAnyJobRunning()); };
  
//   params.appWindowParams.windowTitle = app_name;
//   params.appWindowParams.windowGeometry.windowSizeState = HelloImGui::WindowSizeState::Maximized;
  
//   params.imGuiWindowParams.defaultImGuiWindowType = HelloImGui::DefaultImGuiWindowType::ProvideFullScreenDockSpace;
//   params.imGuiWindowParams.enableViewports = true;
//   params.imGuiWindowParams.showMenuBar = true;
//   params.imGuiWindowParams.showMenu_App = false;
//   params.imGuiWindowParams.showMenu_View = false;
//   params.imGuiWindowParams.showStatusBar = false;  // Use our custom status bar instead
  
//   params.fpsIdling.enableIdling = false;
  
//   // Reserve space for single-line status bar (one frame height + small padding)
//   // The margin is in units relative to frame height, so ~1.0 should be enough for a single line
//   // params.imGuiWindowParams.fullScreenWindow_MarginBottomRight = ImVec2(0.0f, 1.4f);
  
//   params.dockingParams = adsim::CreateDefaultLayout();
    
//   HelloImGui::SetAssetsFolder("assets/");

//   try {
//     HelloImGui::Run(params);
//   } catch(const std::exception& e) {
//     std::cerr << e.what() << '\n';
//   }
  
//   ImPlot::DestroyContext();

//   return 0;
// }
