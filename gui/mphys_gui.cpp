/**
 * mphys_gui.cpp  —  COMSOL-inspired desktop front-end for mphys
 *
 * Layout:
 *   left  (~18%)  Model Builder tree
 *   right (~82%)  Configuration / Results panel
 *
 * Physics modules (1D):
 *   Convection-Diffusion-Reaction  (transient, IDA)
 *   Steady Diffusion               (KINSOL)
 */

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

#include "mphys/boundary_condition.hpp"
#include "mphys/fvm_operators.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/sim_result.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/steady_solver.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/transient_solver.hpp"

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================
// Physics model definitions
// ============================================================

// 1D transient convection-diffusion-reaction:  dc/dt = -u·∂c/∂x + D·∂²c/∂x² - k·c
// BCs: c(0) = c_inlet (Dirichlet), ∂c/∂x(L) = 0 (Neumann)
class ConvDiffModel : public mphys::Model {
 public:
  ConvDiffModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv,
                double D, double u, double k, double c_inlet)
      : Model(mesh, sv), D_(D), u_(u), k_(k) {
    c_ = AddField("c", 0.0);
    SetBcs(c_, {mphys::DirichletBc(c_inlet), mphys::NeumannBc(0.0)});
  }

  void Residual(double,
                const std::vector<mphys::Field>& y,
                const std::vector<mphys::Field>& ydot,
                const std::vector<double>&,
                const std::vector<double>&,
                std::vector<mphys::Field>& rr,
                std::vector<double>&) override {
    rr[c_] = mphys::fvm::Ddt(ydot[c_])
           + mphys::fvm::Convection(y[c_], u_, mesh_, bcs_[c_])
           - mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_])
           + y[c_] * k_;
  }

 private:
  int c_ = 0;
  double D_, u_, k_;
};

// 1D steady diffusion:  -D·∂²c/∂x² = 0
// BCs: c(0) = c_left, c(L) = c_right (both Dirichlet)
class SteadyDiffModel : public mphys::Model {
 public:
  SteadyDiffModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv,
                  double D, double c_left, double c_right)
      : Model(mesh, sv), D_(D) {
    c_ = AddField("c", 0.5 * (c_left + c_right));
    SetBcs(c_, {mphys::DirichletBc(c_left), mphys::DirichletBc(c_right)});
  }

  void Residual(double,
                const std::vector<mphys::Field>& y,
                const std::vector<mphys::Field>&,
                const std::vector<double>&,
                const std::vector<double>&,
                std::vector<mphys::Field>& rr,
                std::vector<double>&) override {
    rr[c_] = -mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_]);
  }

 private:
  int c_ = 0;
  double D_;
};

// ============================================================
// Application state
// ============================================================

enum class NavNode { Geometry, Physics, Mesh, Study, Results };
enum class PhysicsModule { ConvectionDiffusion, SteadyDiffusion };

struct AppState {
  NavNode       nav    = NavNode::Geometry;
  PhysicsModule physics = PhysicsModule::ConvectionDiffusion;

  // Mesh
  int   n_cells = 100;
  float x_start = 0.0f;
  float x_end   = 1.0f;

  // Convection-diffusion parameters
  float D_cd    = 1e-4f;
  float u_cd    = 0.01f;
  float k_cd    = 0.1f;
  float c_inlet = 1.0f;

  // Steady diffusion parameters
  float D_sd    = 1.0f;
  float c_left  = 0.0f;
  float c_right = 1.0f;

  // Solver / time settings
  float t_end      = 50.0f;
  float dt_initial = 1e-3f;
  float dt_max     = 1.0f;
  float rel_tol    = 1e-6f;
  float abs_tol    = 1e-8f;

  // Results
  mphys::SimResult     result;
  std::vector<double>  cell_centres;
  int                  snapshot_idx = 0;
  bool                 has_results  = false;
  std::string          status_msg;
};

// ============================================================
// Simulation runner (synchronous, called from GUI)
// ============================================================

static void RunSimulation(AppState& s) {
  s.has_results = false;
  s.result.snapshots.clear();
  s.cell_centres.clear();
  s.status_msg = "Running...";

  try {
    if (s.x_end <= s.x_start)
      throw std::runtime_error("x_end must be greater than x_start");

    auto mesh = mphys::MakeUniformMesh1D(
        static_cast<double>(s.x_start),
        static_cast<double>(s.x_end),
        s.n_cells);

    s.cell_centres = mesh.cell_centres;

    mphys::SolverOptions opts;
    opts.tolerance.relative = static_cast<double>(s.rel_tol);
    opts.tolerance.absolute = static_cast<double>(s.abs_tol);
    opts.initial_time_step  = static_cast<double>(s.dt_initial);
    opts.maximum_time_step  = static_cast<double>(s.dt_max);

    mphys::SunContext  sunctx;
    mphys::StateVector sv(mesh.n_cells);

    if (s.physics == PhysicsModule::ConvectionDiffusion) {
      ConvDiffModel model(mesh, sv,
          static_cast<double>(s.D_cd),
          static_cast<double>(s.u_cd),
          static_cast<double>(s.k_cd),
          static_cast<double>(s.c_inlet));

      mphys::TransientSolver solver(model, opts, sunctx);
      solver.Solve(0.0, static_cast<double>(s.t_end),
          [&](double t, const std::vector<mphys::Field>& f,
              const std::vector<double>& a) {
            s.result.Record(t, f, a);
          });
    } else {
      SteadyDiffModel model(mesh, sv,
          static_cast<double>(s.D_sd),
          static_cast<double>(s.c_left),
          static_cast<double>(s.c_right));

      mphys::SteadySolver solver(model, opts, sunctx);
      solver.Solve();
      s.result.Record(0.0, model.fields(), model.algebraics());
    }

    s.snapshot_idx = static_cast<int>(s.result.snapshots.size()) - 1;
    s.has_results  = true;
    s.nav          = NavNode::Results;
    s.status_msg   = "Done — " + std::to_string(s.result.snapshots.size()) + " snapshots";

  } catch (const std::exception& e) {
    s.status_msg = std::string("Error: ") + e.what();
  }
}

// ============================================================
// GUI panels
// ============================================================

static void ShowGeometryPanel() {
  ImGui::SeparatorText("Coordinate System");
  ImGui::Spacing();

  static const char* kSystems[] = {"1D", "2D", "2D Axisymmetric", "3D"};
  for (int i = 0; i < 4; ++i) {
    int dummy = 0;
    if (i != 0) ImGui::BeginDisabled();
    ImGui::RadioButton(kSystems[i], &dummy, 0);
    if (i != 0) ImGui::EndDisabled();
    if (i < 3) ImGui::SameLine();
  }
  ImGui::Spacing();
  ImGui::TextDisabled("2D, 2D Axisymmetric, and 3D — coming soon.");
}

static void ShowPhysicsPanel(AppState& s) {
  ImGui::SeparatorText("Physics Module");
  ImGui::Spacing();

  static const char* kModules[] = {"Convection-Diffusion-Reaction", "Steady Diffusion"};
  int sel = static_cast<int>(s.physics);
  if (ImGui::Combo("Module##phys", &sel, kModules, 2))
    s.physics = static_cast<PhysicsModule>(sel);

  ImGui::Spacing();

  if (s.physics == PhysicsModule::ConvectionDiffusion) {
    ImGui::SeparatorText("Parameters");
    ImGui::InputFloat("Diffusivity D  [m²/s]",  &s.D_cd,    0, 0, "%.2e");
    ImGui::InputFloat("Velocity u  [m/s]",       &s.u_cd,    0, 0, "%.4f");
    ImGui::InputFloat("Reaction rate k  [1/s]",  &s.k_cd,    0, 0, "%.4f");
    ImGui::InputFloat("Inlet concentration c₀",  &s.c_inlet, 0, 0, "%.4f");
    ImGui::Spacing();
    ImGui::TextDisabled("PDE:  dc/dt = -u dc/dx + D d²c/dx² - k c");
    ImGui::TextDisabled("BCs:  c(0) = c_inlet  (Dirichlet)");
    ImGui::TextDisabled("      dc/dx(L) = 0    (Neumann)");
    ImGui::TextDisabled("IC:   c(x,0) = 0");
  } else {
    ImGui::SeparatorText("Parameters");
    ImGui::InputFloat("Diffusivity D  [m²/s]",  &s.D_sd,    0, 0, "%.4f");
    ImGui::InputFloat("Left BC value",           &s.c_left,  0, 0, "%.4f");
    ImGui::InputFloat("Right BC value",          &s.c_right, 0, 0, "%.4f");
    ImGui::Spacing();
    ImGui::TextDisabled("PDE:  -D d²c/dx² = 0  (steady)");
    ImGui::TextDisabled("BCs:  c(0) = c_left, c(L) = c_right  (Dirichlet)");
  }
}

static void ShowMeshPanel(AppState& s) {
  ImGui::SeparatorText("Domain");
  ImGui::InputFloat("x start  [m]", &s.x_start, 0, 0, "%.3f");
  ImGui::InputFloat("x end    [m]", &s.x_end,   0, 0, "%.3f");
  if (s.x_end <= s.x_start) s.x_end = s.x_start + 0.01f;

  ImGui::Spacing();
  ImGui::SeparatorText("Discretisation");
  ImGui::SliderInt("Number of cells", &s.n_cells, 10, 500);

  float dx = (s.x_end - s.x_start) / static_cast<float>(s.n_cells);
  ImGui::Spacing();
  ImGui::Text("Cell size: %.5f m", static_cast<double>(dx));
}

static void ShowStudyPanel(AppState& s) {
  bool transient = (s.physics == PhysicsModule::ConvectionDiffusion);

  ImGui::SeparatorText("Solver");
  ImGui::TextUnformatted(transient ? "Transient  (SUNDIALS IDA)" : "Steady-state  (SUNDIALS KINSOL)");

  if (transient) {
    ImGui::Spacing();
    ImGui::SeparatorText("Time Settings");
    ImGui::InputFloat("End time  [s]",     &s.t_end,      0, 0, "%.3f");
    ImGui::InputFloat("Initial dt  [s]",   &s.dt_initial, 0, 0, "%.2e");
    ImGui::InputFloat("Max dt  [s]",       &s.dt_max,     0, 0, "%.3f");
  }

  ImGui::Spacing();
  ImGui::SeparatorText("Tolerances");
  ImGui::InputFloat("Relative tolerance", &s.rel_tol, 0, 0, "%.2e");
  ImGui::InputFloat("Absolute tolerance", &s.abs_tol, 0, 0, "%.2e");

  ImGui::Spacing();
  ImGui::Spacing();

  bool can_run = !s.status_msg.starts_with("Running");
  if (!can_run) ImGui::BeginDisabled();
  if (ImGui::Button("Run", ImVec2(120, 32)))
    RunSimulation(s);
  if (!can_run) ImGui::EndDisabled();

  if (!s.status_msg.empty()) {
    ImGui::SameLine();
    ImGui::TextUnformatted(s.status_msg.c_str());
  }
}

static void ShowResultsPanel(AppState& s) {
  if (!s.has_results || s.result.snapshots.empty()) {
    ImGui::TextDisabled("No results yet — run a simulation from the Study node.");
    return;
  }

  const auto& snap = s.result.snapshots[s.snapshot_idx];
  int n_snaps = static_cast<int>(s.result.snapshots.size());

  // Time slider (transient only)
  if (n_snaps > 1) {
    ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - 120.0f);
    ImGui::SliderInt("Snapshot", &s.snapshot_idx, 0, n_snaps - 1);
    ImGui::SameLine();
    ImGui::Text("t = %.4f s", snap.t);
  } else {
    ImGui::Text("Steady-state result");
  }

  ImGui::Spacing();

  // Field selector
  int n_fields = static_cast<int>(snap.fields.size());
  static int field_idx = 0;
  if (field_idx >= n_fields) field_idx = 0;

  if (n_fields > 1) {
    std::vector<const char*> names;
    names.reserve(n_fields);
    for (const auto& f : snap.fields) names.push_back(f.name.c_str());
    ImGui::Combo("Field", &field_idx, names.data(), n_fields);
    ImGui::Spacing();
  }

  // Plot
  const auto& field = snap.fields[field_idx];
  int n_cells = static_cast<int>(s.cell_centres.size());

  ImVec2 plot_size = ImGui::GetContentRegionAvail();
  plot_size.y -= 4.0f;

  if (ImPlot::BeginPlot("##results", plot_size)) {
    ImPlot::SetupAxes("x  [m]", field.name.c_str());
    ImPlot::SetupAxisLimitsConstraints(ImAxis_X1,
        static_cast<double>(s.x_start),
        static_cast<double>(s.x_end));

    ImPlot::PlotLine(field.name.c_str(),
        s.cell_centres.data(),
        field.values.data(),
        n_cells);

    ImPlot::EndPlot();
  }
}

static void ShowConfigPanel(AppState& s) {
  switch (s.nav) {
    case NavNode::Geometry: ShowGeometryPanel();   break;
    case NavNode::Physics:  ShowPhysicsPanel(s);   break;
    case NavNode::Mesh:     ShowMeshPanel(s);       break;
    case NavNode::Study:    ShowStudyPanel(s);      break;
    case NavNode::Results:  ShowResultsPanel(s);    break;
  }
}

// ============================================================
// Dock layout (set up once on first frame)
// ============================================================

static void SetupDockLayout(ImGuiID dockspace_id) {
  ImGui::DockBuilderRemoveNode(dockspace_id);
  ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
  ImGui::DockBuilderSetNodeSize(dockspace_id, ImGui::GetMainViewport()->WorkSize);

  ImGuiID left, right;
  ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.18f, &left, &right);

  ImGui::DockBuilderDockWindow("Model Builder", left);
  ImGui::DockBuilderDockWindow("Configuration", right);
  ImGui::DockBuilderFinish(dockspace_id);
}

// ============================================================
// Main render function — called each frame
// ============================================================

static void RenderFrame(AppState& s) {
  static bool first_frame = true;

  // Full-screen DockSpace
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

  // Menu bar
  if (ImGui::BeginMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem("Quit", "Alt+F4"))
        glfwSetWindowShouldClose(glfwGetCurrentContext(), GLFW_TRUE);
      ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Study")) {
      if (ImGui::MenuItem("Run", "F5"))
        RunSimulation(s);
      ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Help")) {
      if (ImGui::BeginMenu("About")) {
        ImGui::TextUnformatted("mphys GUI — COMSOL-inspired 1D physics front-end");
        ImGui::TextUnformatted("Powered by Dear ImGui + ImPlot + SUNDIALS");
        ImGui::EndMenu();
      }
      ImGui::EndMenu();
    }

    // Quick-run button in menu bar
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + 12.0f);
    if (ImGui::SmallButton("  Run  "))
      RunSimulation(s);

    if (!s.status_msg.empty()) {
      ImGui::SetCursorPosX(ImGui::GetCursorPosX() + 8.0f);
      ImGui::TextDisabled("%s", s.status_msg.c_str());
    }

    ImGui::EndMenuBar();
  }

  // DockSpace
  ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
  ImGui::DockSpace(dockspace_id, ImVec2(0, 0), ImGuiDockNodeFlags_None);

  if (first_frame) {
    first_frame = false;
    SetupDockLayout(dockspace_id);
  }

  ImGui::End();

  // ---- Model Builder window ----
  ImGuiWindowFlags panel_flags = ImGuiWindowFlags_NoTitleBar;
  ImGui::Begin("Model Builder", nullptr, 0);
  ImGui::TextColored(ImVec4(0.4f, 0.6f, 1.0f, 1.0f), "Model Builder");
  ImGui::Separator();
  ImGui::Spacing();

  static const char* kNodes[] = {"Geometry", "Physics", "Mesh", "Study", "Results"};
  for (int i = 0; i < 5; ++i) {
    NavNode node = static_cast<NavNode>(i);
    bool selected = (s.nav == node);

    ImGuiTreeNodeFlags flags =
        ImGuiTreeNodeFlags_Leaf |
        ImGuiTreeNodeFlags_NoTreePushOnOpen |
        ImGuiTreeNodeFlags_SpanFullWidth;
    if (selected) flags |= ImGuiTreeNodeFlags_Selected;

    if (i == 4 && s.has_results) {
      ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.2f, 0.75f, 0.3f, 1.0f));
    }
    ImGui::TreeNodeEx(kNodes[i], flags);
    if (i == 4 && s.has_results) ImGui::PopStyleColor();

    if (ImGui::IsItemClicked()) s.nav = node;
  }

  ImGui::End();  // Model Builder

  // ---- Configuration / Results window ----
  ImGui::Begin("Configuration", nullptr, 0);

  // Breadcrumb header
  static const char* kNodeTitles[] = {
    "Geometry", "Physics", "Mesh", "Study", "Results"
  };
  ImGui::TextColored(ImVec4(0.4f, 0.6f, 1.0f, 1.0f),
      "%s", kNodeTitles[static_cast<int>(s.nav)]);
  ImGui::Separator();
  ImGui::Spacing();

  ShowConfigPanel(s);

  ImGui::End();  // Configuration
}

// ============================================================
// Entry point
// ============================================================

int main() {
  glfwSetErrorCallback([](int error, const char* desc) {
    fprintf(stderr, "GLFW error %d: %s\n", error, desc);
  });

  if (!glfwInit()) return 1;

#ifdef __APPLE__
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

  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

  ImGui::StyleColorsDark();
  ImGui::GetStyle().FrameRounding    = 3.0f;
  ImGui::GetStyle().WindowRounding   = 4.0f;
  ImGui::GetStyle().GrabRounding     = 3.0f;
  ImGui::GetStyle().TabRounding      = 3.0f;
  ImGui::GetStyle().WindowBorderSize = 1.0f;

  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  AppState state;

  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    RenderFrame(state);

    ImGui::Render();

    int fb_w, fb_h;
    glfwGetFramebufferSize(window, &fb_w, &fb_h);
    glViewport(0, 0, fb_w, fb_h);
    glClearColor(0.12f, 0.12f, 0.12f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(window);
  }

  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImPlot::DestroyContext();
  ImGui::DestroyContext();

  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}
