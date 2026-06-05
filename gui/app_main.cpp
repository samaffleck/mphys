/**
 * app_main.cpp  —  COMSOL-inspired front-end for mphys (entry point + shell)
 *
 * Layout:
 *   left  (~18%)  Model Builder tree
 *   centre (~42%) Configuration panel
 *   right         Geometry View / Results canvas
 *
 * The application is split across several translation units:
 *   app_state.hpp     plain data structures (no UI)
 *   expr.hpp          arithmetic expression evaluator
 *   gui_widgets.*     reusable ImGui input widgets
 *   simulation.*      mesh building, model config, solver runners
 *   serialization.*   JSON save/load
 *   panels_config.*   Geometry / Physics / Mesh / Study forms
 *   panels_view.*     Geometry preview + result visualisation
 *   app_main.cpp      window/loop, menus, dock layout (this file)
 */

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include <GLFW/glfw3.h>
#include "themes.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/html5.h>
#else
#include "tinyfiledialogs.h"
#endif

#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"
#include "implot3d.h"

#include "app_state.hpp"
#include "panels.hpp"
#include "serialization.hpp"
#include "simulation.hpp"

// ============================================================
// Config panel dispatcher
// ============================================================

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
        static const std::string kExDir = std::string(MPHYS_ASSETS_DIR) + "/examples/";
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
        ImGui::TextUnformatted("mphys GUI — COMSOL-inspired 1D/2D/3D physics front-end");
        ImGui::TextUnformatted("Powered by Dear ImGui + ImPlot + ImPlot3D + SUNDIALS");
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

  GLFWwindow* window = glfwCreateWindow(1400, 900, "mphys", nullptr, nullptr);
  if (!window) { glfwTerminate(); return 1; }

  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImPlot::CreateContext();
  ImPlot3D::CreateContext();

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
  ImGui_ImplGlfw_InstallEmscriptenCallbacks(window, "#canvas");
#endif
  ImGui_ImplOpenGL3_Init(glsl_version);

#if defined(__EMSCRIPTEN__)
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
  ImPlot3D::DestroyContext();
  ImPlot::DestroyContext();
  ImGui::DestroyContext();
  glfwDestroyWindow(window);
  glfwTerminate();
#endif
  return 0;
}
