#include "panels.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "imgui.h"
#include "implot.h"
#include "implot3d.h"

#include "gui_widgets.hpp"
#include "mphys/mesh.hpp"

// ============================================================
// Geometry View — CAD canvas
// ============================================================

// Structured-box preview for the 2D/3D convection-diffusion-reaction model.
// Boundary edges/faces are tinted by BC type (Dirichlet = orange, Neumann =
// blue) so the inflow/outflow setup is visible at a glance.
static void ShowBoxPreview(const AppState& s, ImDrawList* dl, ImVec2 orig,
                           float cw, float ch) {
  const BoxGeometry& b = s.box;
  auto bc_col = [&](int p) {
    return b.bc_type[p] == 0 ? IM_COL32(255, 170, 90, 255)
                             : IM_COL32(90, 160, 255, 255);
  };

  char title[80];
  if (b.dim == 2) {
    snprintf(title, sizeof(title), "2D domain   %d x %d cells", b.nx, b.ny);

    float Lx = b.x1 - b.x0, Ly = b.y1 - b.y0;
    if (Lx <= 0.0f) Lx = 1.0f;
    if (Ly <= 0.0f) Ly = 1.0f;
    const float pad = 80.0f;
    float sc = std::min((cw - 2 * pad) / Lx, (ch - 2 * pad) / Ly);
    float rw = Lx * sc, rh = Ly * sc;
    ImVec2 p0(orig.x + (cw - rw) * 0.5f, orig.y + (ch - rh) * 0.5f);
    ImVec2 p1(p0.x + rw, p0.y + rh);

    dl->AddRectFilled(p0, p1, IM_COL32(70, 100, 150, 40), 0.0f);
    // Edges. Screen top = y max, screen bottom = y min.
    dl->AddLine(ImVec2(p0.x, p0.y), ImVec2(p0.x, p1.y), bc_col(0), 3.0f);  // left
    dl->AddLine(ImVec2(p1.x, p0.y), ImVec2(p1.x, p1.y), bc_col(1), 3.0f);  // right
    dl->AddLine(ImVec2(p0.x, p1.y), ImVec2(p1.x, p1.y), bc_col(2), 3.0f);  // bottom
    dl->AddLine(ImVec2(p0.x, p0.y), ImVec2(p1.x, p0.y), bc_col(3), 3.0f);  // top

    char lx[32]; snprintf(lx, sizeof(lx), "x: %.3g m", (double)Lx);
    char ly[32]; snprintf(ly, sizeof(ly), "y: %.3g m", (double)Ly);
    dl->AddText(ImVec2((p0.x + p1.x) * 0.5f - 30, p1.y + 8),
                IM_COL32(170, 185, 210, 220), lx);
    dl->AddText(ImVec2(p0.x - 56, (p0.y + p1.y) * 0.5f - 8),
                IM_COL32(170, 185, 210, 220), ly);
  } else {
    snprintf(title, sizeof(title), "3D domain   %d x %d x %d cells", b.nx, b.ny,
             b.nz);

    float Lx = b.x1 - b.x0, Ly = b.y1 - b.y0, Lz = b.z1 - b.z0;
    float mx = std::max({Lx, Ly, Lz, 1e-6f});
    float ax = Lx / mx, ay = Ly / mx, az = Lz / mx;
    float sc = std::min(cw, ch) * 0.34f;
    ImVec2 oc(orig.x + cw * 0.5f, orig.y + ch * 0.55f);
    auto proj = [&](float x, float y, float z) -> ImVec2 {
      x -= ax * 0.5f; y -= ay * 0.5f; z -= az * 0.5f;
      return ImVec2(oc.x + sc * (x * 0.866f - z * 0.866f),
                    oc.y + sc * (-y + x * 0.5f + z * 0.5f));
    };
    ImVec2 c[8];
    for (int i = 0; i < 8; ++i)
      c[i] = proj((i & 1) ? ax : 0.0f, (i & 2) ? ay : 0.0f, (i & 4) ? az : 0.0f);
    static const int edges[12][2] = {
        {0,1},{2,3},{4,5},{6,7}, {0,2},{1,3},{4,6},{5,7}, {0,4},{1,5},{2,6},{3,7}};
    for (auto& e : edges)
      dl->AddLine(c[e[0]], c[e[1]], IM_COL32(120, 160, 230, 220), 2.0f);
    for (int i = 0; i < 8; ++i)
      dl->AddCircleFilled(c[i], 2.5f, IM_COL32(200, 220, 255, 230));

    dl->AddText(ImVec2(proj(ax * 0.5f, 0, 0).x - 8, proj(ax * 0.5f, 0, 0).y + 6),
                IM_COL32(170, 185, 210, 220), "x");
    dl->AddText(ImVec2(proj(0, ay * 0.5f, 0).x - 16, proj(0, ay * 0.5f, 0).y),
                IM_COL32(170, 185, 210, 220), "y");
    dl->AddText(ImVec2(proj(0, 0, az * 0.5f).x - 4, proj(0, 0, az * 0.5f).y + 6),
                IM_COL32(170, 185, 210, 220), "z");
  }

  ImVec2 tts = ImGui::CalcTextSize(title);
  dl->AddText(ImVec2(orig.x + (cw - tts.x) * 0.5f, orig.y + 16.0f),
              IM_COL32(150, 170, 200, 220), title);
}

void ShowGeometryView(AppState& s) {
  static constexpr const char* kUnits[] = {"m", "cm", "mm", "\xc2\xb5m"};
  Geometry1D& geo = s.geo;

  ImDrawList* dl   = ImGui::GetWindowDrawList();
  ImVec2      orig = ImGui::GetWindowPos();
  float       cw   = ImGui::GetWindowWidth();
  float       ch   = ImGui::GetWindowHeight();

  // 2D/3D structured-box preview.
  if (s.model_id == kConvDiffId && s.dim >= 2) {
    ShowBoxPreview(s, dl, orig, cw, ch);
    return;
  }

  if (s.model_id == kSpmId) {
    // Schematic of the two representative particles, sized by relative radius.
    float cy   = orig.y + ch * 0.5f;
    float rmax = std::max(s.spm.R_n, s.spm.R_p);
    if (rmax <= 0.0f) rmax = 1.0f;
    float base = std::min(cw, ch) * 0.18f;
    float rn   = base * (s.spm.R_n / rmax);
    float rp   = base * (s.spm.R_p / rmax);
    float cxn  = orig.x + cw * 0.33f;
    float cxp  = orig.x + cw * 0.67f;

    auto draw_particle = [&](float cx, float r, ImU32 col, const char* tag,
                             float radius_m) {
      for (int k = 4; k >= 1; --k) {  // concentric shells → diffusion hint
        float rr = r * k / 4.0f;
        dl->AddCircle(ImVec2(cx, cy), rr,
                      IM_COL32(((col >> 0) & 0xFF), ((col >> 8) & 0xFF),
                               ((col >> 16) & 0xFF), 60), 0, 1.0f);
      }
      dl->AddCircle(ImVec2(cx, cy), r, col, 0, 2.5f);
      dl->AddCircleFilled(ImVec2(cx, cy), 2.5f, col);
      char lbl[64];
      snprintf(lbl, sizeof(lbl), "%s   R = %.3g m", tag, (double)radius_m);
      ImVec2 ts = ImGui::CalcTextSize(lbl);
      dl->AddText(ImVec2(cx - ts.x * 0.5f, cy + r + 10.0f),
                  IM_COL32(190, 205, 230, 220), lbl);
    };

    draw_particle(cxn, rn, IM_COL32(120, 185, 255, 255), "Negative", s.spm.R_n);
    draw_particle(cxp, rp, IM_COL32(255, 170, 90, 255), "Positive", s.spm.R_p);

    const char* title = "Single Particle Model";
    ImVec2 tts = ImGui::CalcTextSize(title);
    dl->AddText(ImVec2(orig.x + (cw - tts.x) * 0.5f, orig.y + 16.0f),
                IM_COL32(150, 170, 200, 200), title);
    return;
  }

  if (s.model_id == kSpmeId) {
    // Cell sandwich (negative | separator | positive) with the two
    // representative particles drawn above their electrodes.
    const SpmeInputs& m = s.spme;
    float Ln = m.core.L_n, Ls = m.L_s, Lp = m.core.L_p;
    float Ltot = std::max(Ln + Ls + Lp, 1e-12f);

    float pad   = 60.0f;
    float usable = cw - 2.0f * pad;
    float x0px  = orig.x + pad;
    float bandTop = orig.y + ch * 0.56f;
    float bandBot = orig.y + ch * 0.78f;
    float xn = x0px + usable * (Ln / Ltot);
    float xs = x0px + usable * ((Ln + Ls) / Ltot);
    float xe = x0px + usable;

    auto band = [&](float xa, float xb, ImU32 col, const char* lbl) {
      dl->AddRectFilled(ImVec2(xa, bandTop), ImVec2(xb, bandBot), col, 3.0f);
      dl->AddRect(ImVec2(xa, bandTop), ImVec2(xb, bandBot),
                  IM_COL32(200, 210, 230, 120), 3.0f);
      ImVec2 ts = ImGui::CalcTextSize(lbl);
      if (xb - xa > ts.x + 6.0f)
        dl->AddText(ImVec2((xa + xb) * 0.5f - ts.x * 0.5f,
                           (bandTop + bandBot) * 0.5f - ts.y * 0.5f),
                    IM_COL32(20, 25, 35, 230), lbl);
    };
    band(x0px, xn, IM_COL32(120, 185, 255, 200), "Negative");
    band(xn,   xs, IM_COL32(170, 180, 200, 180), "Separator");
    band(xs,   xe, IM_COL32(255, 170, 90, 200),  "Positive");

    // Particles above each electrode.
    float rmax = std::max(m.core.R_n, m.core.R_p);
    if (rmax <= 0.0f) rmax = 1.0f;
    float base = std::min(cw, ch) * 0.11f;
    float pcy  = orig.y + ch * 0.30f;
    auto particle = [&](float cx, float radius_m, ImU32 col) {
      float r = base * (radius_m / rmax);
      for (int k = 3; k >= 1; --k)
        dl->AddCircle(ImVec2(cx, pcy), r * k / 3.0f,
                      IM_COL32(((col >> 0) & 0xFF), ((col >> 8) & 0xFF),
                               ((col >> 16) & 0xFF), 60), 0, 1.0f);
      dl->AddCircle(ImVec2(cx, pcy), r, col, 0, 2.5f);
      dl->AddCircleFilled(ImVec2(cx, pcy), 2.0f, col);
    };
    particle((x0px + xn) * 0.5f, m.core.R_n, IM_COL32(120, 185, 255, 255));
    particle((xs + xe) * 0.5f,   m.core.R_p, IM_COL32(255, 170, 90, 255));

    dl->AddText(ImVec2(x0px, bandBot + 8.0f), IM_COL32(150, 170, 200, 200),
                "x = 0");
    const char* xl = "x = L";
    ImVec2 xls = ImGui::CalcTextSize(xl);
    dl->AddText(ImVec2(xe - xls.x, bandBot + 8.0f),
                IM_COL32(150, 170, 200, 200), xl);

    const char* title = "Single Particle Model w/ Electrolyte";
    ImVec2 tts = ImGui::CalcTextSize(title);
    dl->AddText(ImVec2(orig.x + (cw - tts.x) * 0.5f, orig.y + 16.0f),
                IM_COL32(150, 170, 200, 200), title);
    return;
  }

  if (!geo.built || (int)geo.nodes.size() < 2) {
    const char* msg = "Configure geometry and press  Build";
    ImVec2 ts = ImGui::CalcTextSize(msg);
    dl->AddText(ImVec2(orig.x + (cw - ts.x) * 0.5f, orig.y + (ch - ts.y) * 0.5f),
                IM_COL32(100, 100, 120, 200), msg);
    return;
  }

  const float pad    = 60.0f;
  const float usable = cw - 2.0f * pad;
  const float line_y = orig.y + ch * 0.50f;

  float x_min = geo.nodes.front().x;
  float x_max = geo.nodes.back().x;
  float span  = x_max - x_min;
  if (span <= 0.0f) span = 1.0f;

  auto node_px = [&](float x) -> float {
    return orig.x + pad + (x - x_min) / span * usable;
  };

  ImVec2 mouse     = ImGui::GetMousePos();
  int hovered_node   = -1;
  int hovered_domain = -1;
  bool canvas_hov  = ImGui::IsWindowHovered(ImGuiHoveredFlags_ChildWindows);

  if (canvas_hov) {
    for (int i = 0; i < (int)geo.nodes.size(); ++i) {
      float px = node_px(geo.nodes[i].x);
      float dx = mouse.x - px, dy = mouse.y - line_y;
      if (dx * dx + dy * dy < 9.0f * 9.0f) { hovered_node = i; break; }
    }
    if (hovered_node < 0) {
      for (int d = 0; d < (int)geo.domains.size(); ++d) {
        float px_l = node_px(geo.nodes[d].x);
        float px_r = node_px(geo.nodes[d + 1].x);
        if (mouse.x > px_l && mouse.x < px_r && fabsf(mouse.y - line_y) < 7.0f)
          hovered_domain = d;
      }
    }
  }

  if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
    if (hovered_node >= 0) {
      geo.selected_node = hovered_node; geo.selected_domain = -1;
    } else if (hovered_domain >= 0) {
      geo.selected_domain = hovered_domain; geo.selected_node = -1;
    } else if (canvas_hov) {
      geo.selected_node = geo.selected_domain = -1;
    }
  }

  if (hovered_node >= 0 || hovered_domain >= 0)
    ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);

  // Faint axis
  dl->AddLine(ImVec2(orig.x + pad * 0.5f, line_y),
              ImVec2(orig.x + cw - pad * 0.5f, line_y),
              IM_COL32(80, 80, 100, 60), 1.0f);

  // Domain segments
  for (int d = 0; d < (int)geo.domains.size(); ++d) {
    float px_l = node_px(geo.nodes[d].x);
    float px_r = node_px(geo.nodes[d + 1].x);
    bool  sel  = (d == geo.selected_domain);
    bool  hov  = (d == hovered_domain);
    ImU32 col  = sel ? IM_COL32(255, 190, 50, 255) :
                 hov ? IM_COL32(120, 185, 255, 255) :
                       IM_COL32(75,  145, 245, 255);
    dl->AddLine(ImVec2(px_l, line_y), ImVec2(px_r, line_y), col, 3.5f);

    char lbl[8]; snprintf(lbl, sizeof(lbl), "%d", d + 1);
    float mid_x = (px_l + px_r) * 0.5f;
    ImVec2 ts = ImGui::CalcTextSize(lbl);
    dl->AddText(ImVec2(mid_x - ts.x * 0.5f, line_y - ts.y - 10.0f),
                sel ? IM_COL32(255, 210, 100, 200) : IM_COL32(120, 160, 220, 160), lbl);
  }

  // Nodes
  for (int i = 0; i < (int)geo.nodes.size(); ++i) {
    float px       = node_px(geo.nodes[i].x);
    bool  sel      = (i == geo.selected_node);
    bool  hov      = (i == hovered_node);
    bool  terminal = (i == 0 || i == (int)geo.nodes.size() - 1);
    float r        = terminal ? 6.5f : 4.5f;

    ImU32 fill = sel ? IM_COL32(255, 200, 60, 255) :
                 hov ? IM_COL32(220, 235, 255, 255) :
                       IM_COL32(200, 220, 255, terminal ? 255 : 180);
    dl->AddCircleFilled(ImVec2(px, line_y), r, fill);
    if (sel)
      dl->AddCircle(ImVec2(px, line_y), r + 3.5f, IM_COL32(255, 200, 60, 180), 0, 1.5f);
    else if (hov)
      dl->AddCircle(ImVec2(px, line_y), r + 2.5f, IM_COL32(200, 220, 255, 120), 0, 1.0f);

    dl->AddLine(ImVec2(px, line_y + r + 1.0f), ImVec2(px, line_y + r + 5.0f),
                IM_COL32(150, 170, 200, 140), 1.0f);

    if (terminal || sel || hov) {
      char xlbl[48]; snprintf(xlbl, sizeof(xlbl), "%.4g %s",
                               (double)geo.nodes[i].x, kUnits[s.geo_unit]);
      ImVec2 xs = ImGui::CalcTextSize(xlbl);
      ImU32  tc = sel ? IM_COL32(255, 210, 100, 220) : IM_COL32(170, 185, 210, 200);
      dl->AddText(ImVec2(px - xs.x * 0.5f, line_y + r + 7.0f), tc, xlbl);
    }
  }

  // Selection overlay
  if (geo.selected_node >= 0 || geo.selected_domain >= 0) {
    const char* u = kUnits[s.geo_unit];
    char info[128] = {};
    if (geo.selected_node >= 0) {
      int  i    = geo.selected_node;
      bool term = (i == 0 || i == (int)geo.nodes.size() - 1);
      snprintf(info, sizeof(info), "Node %d  x = %.6g %s%s",
               i, (double)geo.nodes[i].x, u, term ? "" : "  (internal)");
    } else {
      int d = geo.selected_domain;
      double L = geo.nodes[d + 1].x - geo.nodes[d].x;
      snprintf(info, sizeof(info), "Domain %d  L = %.6g %s  n = %d",
               d + 1, L, u, geo.domains[d].n_cells);
    }
    ImVec2 ts = ImGui::CalcTextSize(info);
    float  ix = orig.x + 10.0f, iy = orig.y + ch - ts.y - 10.0f;
    dl->AddRectFilled(ImVec2(ix - 4, iy - 3), ImVec2(ix + ts.x + 4, iy + ts.y + 3),
                      IM_COL32(20, 20, 30, 180), 3.0f);
    dl->AddText(ImVec2(ix, iy), IM_COL32(200, 210, 230, 230), info);
  }
}

// ============================================================
// Results panel
// ============================================================

// Linearly interpolate a 1D field at time t between stored snapshots.
static std::vector<double> InterpolateField(const mphys::SimResult& res,
                                             int field_idx, double t) {
  const auto& snaps = res.snapshots;
  if (snaps.empty()) return {};
  if (t <= snaps.front().t) return snaps.front().fields[field_idx].values;
  if (t >= snaps.back().t)  return snaps.back().fields[field_idx].values;

  int lo = 0, hi = static_cast<int>(snaps.size()) - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (snaps[mid].t <= t) lo = mid; else hi = mid;
  }
  double t0 = snaps[lo].t, t1 = snaps[hi].t;
  double alpha = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0;
  const auto& f0 = snaps[lo].fields[field_idx].values;
  const auto& f1 = snaps[hi].fields[field_idx].values;
  std::vector<double> out(f0.size());
  for (int i = 0; i < (int)f0.size(); ++i)
    out[i] = f0[i] * (1.0 - alpha) + f1[i] * alpha;
  return out;
}

// Linearly interpolate a mesh field (single component) at time t.
static std::vector<double> InterpolateMeshField(
    const std::vector<MeshFieldSnapshot>& snaps, double t) {
  if (snaps.empty()) return {};
  if (t <= snaps.front().t) return snaps.front().values;
  if (t >= snaps.back().t)  return snaps.back().values;
  int lo = 0, hi = static_cast<int>(snaps.size()) - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (snaps[mid].t <= t) lo = mid; else hi = mid;
  }
  double t0 = snaps[lo].t, t1 = snaps[hi].t;
  double a = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0;
  const auto& f0 = snaps[lo].values;
  const auto& f1 = snaps[hi].values;
  std::vector<double> out(f0.size());
  for (size_t i = 0; i < f0.size(); ++i) out[i] = f0[i] * (1.0 - a) + f1[i] * a;
  return out;
}

// Shared time-scrubber for transient runs. Returns true if there is more than
// one snapshot (i.e. the result is transient).
static bool TimeScrubber(AppState& s, double t_min, double t_max, int n_snaps) {
  if (n_snaps <= 1) {
    ImGui::TextUnformatted("Steady-state result");
    ImGui::Spacing();
    return false;
  }
  s.plot_time = std::clamp(s.plot_time, (float)t_min, (float)t_max);
  float avail = ImGui::GetContentRegionAvail().x;
  ImGui::SetNextItemWidth(avail - kUnitColWidth);
  ImGui::SliderFloat("##tslider", &s.plot_time, (float)t_min, (float)t_max, "");
  ImGui::SameLine(0, 6);
  ImGui::TextDisabled("s");
  ImGui::Spacing();
  float w = avail - ImGui::CalcTextSize("t =").x - 6.0f - kUnitColWidth;
  ImGui::TextUnformatted("t ="); ImGui::SameLine(0, 6);
  ImGui::SetNextItemWidth(w);
  ExprInputFloat("##texact", &s.plot_time);
  s.plot_time = std::clamp(s.plot_time, (float)t_min, (float)t_max);
  ImGui::SameLine(0, 6); ImGui::TextDisabled("s");
  ImGui::Spacing();
  ImGui::TextDisabled("%d snapshots  [%.4g s … %.4g s]", n_snaps, t_min, t_max);
  ImGui::Spacing();
  return true;
}

// 2D / 3D field visualisation: heatmap (2D) or coloured point cloud (3D).
static void ShowMeshResults(AppState& s) {
  const auto& snaps = s.mesh_snaps;
  if (snaps.empty()) {
    ImGui::TextDisabled("No results yet — run a simulation from the Study node.");
    return;
  }
  const mphys::Mesh& mesh = s.mesh_result;

  ImGui::SeparatorText(s.result_dim == 3 ? "Field c  (3D)" : "Field c  (2D)");
  ImGui::Spacing();

  const int n_snaps = static_cast<int>(snaps.size());
  const bool transient =
      TimeScrubber(s, snaps.front().t, snaps.back().t, n_snaps);

  std::vector<double> vals = transient
      ? InterpolateMeshField(snaps, static_cast<double>(s.plot_time))
      : snaps.front().values;
  if (vals.empty()) return;

  double vmin = vals[0], vmax = vals[0];
  for (double v : vals) { vmin = std::min(vmin, v); vmax = std::max(vmax, v); }
  if (vmax - vmin < 1e-12) { vmax = vmin + 1.0; }

  const auto& si = mesh.structured;
  const double x0 = si.x0, x1 = si.x0 + si.nx * si.dx;
  const double y0 = si.y0, y1 = si.y0 + si.ny * si.dy;

  const float scale_w = 90.0f;
  ImVec2 avail = ImGui::GetContentRegionAvail();
  ImVec2 psz(std::max(avail.x - scale_w, 50.0f), std::max(avail.y - 6.0f, 50.0f));

  if (s.result_dim == 2) {
    // Reorder row-major (index = j*nx + i) so ImPlot draws y increasing upward
    // (its row 0 maps to the top of the plot, i.e. y max).
    const int nx = si.nx, ny = si.ny;
    std::vector<double> heat(static_cast<size_t>(nx) * ny);
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        heat[static_cast<size_t>(ny - 1 - j) * nx + i] = vals[j * nx + i];

    ImPlot::PushColormap(ImPlotColormap_Viridis);
    if (ImPlot::BeginPlot("##heat", psz, ImPlotFlags_NoLegend | ImPlotFlags_Equal)) {
      ImPlot::SetupAxes("x  [m]", "y  [m]");
      ImPlot::SetupAxesLimits(x0, x1, y0, y1, ImPlotCond_Once);
      ImPlot::PlotHeatmap("c", heat.data(), ny, nx, vmin, vmax, nullptr,
                          ImPlotPoint(x0, y0), ImPlotPoint(x1, y1));
      ImPlot::EndPlot();
    }
    ImGui::SameLine();
    ImPlot::ColormapScale("c", vmin, vmax, ImVec2(scale_w - 16.0f, psz.y));
    ImPlot::PopColormap();
  } else {
    // 3D point cloud: one marker per cell centroid, coloured by value.
    const int nc = mesh.NCells();
    std::vector<float> xs(nc), ys(nc), zs(nc);
    std::vector<ImU32> cols(nc);
    const double inv = 1.0 / (vmax - vmin);
    for (int c = 0; c < nc; ++c) {
      const auto& p = mesh.cells[c].centroid;
      xs[c] = static_cast<float>(p[0]);
      ys[c] = static_cast<float>(p[1]);
      zs[c] = static_cast<float>(p[2]);
      float t = static_cast<float>(std::clamp((vals[c] - vmin) * inv, 0.0, 1.0));
      cols[c] = ImGui::ColorConvertFloat4ToU32(
          ImPlot3D::SampleColormap(t, ImPlot3DColormap_Viridis));
    }
    const double z0 = si.z0, z1 = si.z0 + si.nz * si.dz;

    if (ImPlot3D::BeginPlot("##vol", psz, ImPlot3DFlags_NoLegend)) {
      ImPlot3D::SetupAxes("x", "y", "z");
      ImPlot3D::SetupAxesLimits(x0, x1, y0, y1, z0, z1, ImPlot3DCond_Once);
      ImPlot3DSpec spec;
      spec.Marker = ImPlot3DMarker_Square;
      spec.MarkerSize = 2.5f;
      spec.MarkerFillColors = cols.data();
      spec.MarkerLineColors = cols.data();
      ImPlot3D::PlotScatter("c", xs.data(), ys.data(), zs.data(), nc, spec);
      ImPlot3D::EndPlot();
    }
    ImGui::SameLine();
    ImPlot::ColormapScale("c", vmin, vmax, ImVec2(scale_w - 16.0f, psz.y),
                          "%g", 0, ImPlotColormap_Viridis);
  }
}

void ShowResultsPanel(AppState& s) {
  // 2D/3D mesh results render through their own visualisation.
  if (s.has_results && s.result_dim >= 2) { ShowMeshResults(s); return; }

  if (!s.has_results || s.result.snapshots.empty()) {
    ImGui::TextDisabled("No results yet — run a simulation from the Study node.");
    return;
  }

  const bool is_spm  = (s.model_id == kSpmId);
  const bool is_spme = (s.model_id == kSpmeId);
  const bool is_particle = is_spm || is_spme;

  // SPM/SPMe: terminal voltage vs time, with a marker at the selected instant.
  if (is_particle && !s.spm_voltage.empty()) {
    ImGui::SeparatorText("Terminal Voltage");
    ImGui::Spacing();
    int nv = static_cast<int>(s.spm_voltage.size());
    double v_now = s.spm_voltage.back();
    {
      int best = nv - 1; double bd = 1e300;
      for (int i = 0; i < nv; ++i) {
        double d = std::abs(s.spm_time[i] - (double)s.plot_time);
        if (d < bd) { bd = d; best = i; }
      }
      v_now = s.spm_voltage[best];
    }
    ImGui::Text("V(t = %.4g s) = %.4f V", (double)s.plot_time, v_now);
    ImGui::Spacing();

    ImVec2 vsz = ImGui::GetContentRegionAvail();
    vsz.y *= 0.42f;
    if (ImPlot::BeginPlot("##spm_voltage", vsz)) {
      ImPlot::SetupAxes("t  [s]", "V  [V]",
                        ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
      ImPlot::PlotLine("V", s.spm_time.data(), s.spm_voltage.data(), nv);
      double tx = (double)s.plot_time, vy = v_now;
      ImPlot::PlotScatter("##now", &tx, &vy, 1);
      ImPlot::EndPlot();
    }
    ImGui::Spacing();
    ImGui::SeparatorText(is_spme ? "Concentration Profiles"
                                 : "Particle Concentration");
    ImGui::Spacing();
  }

  const auto& snaps   = s.result.snapshots;
  int         n_snaps = static_cast<int>(snaps.size());
  bool        transient = n_snaps > 1;

  int n_fields = static_cast<int>(snaps[0].fields.size());
  static int field_idx = 0;
  if (field_idx >= n_fields) field_idx = 0;

  if (!transient) {
    ImGui::TextUnformatted("Steady-state result");
    ImGui::Spacing();
  } else {
    float t_min = static_cast<float>(snaps.front().t);
    float t_max = static_cast<float>(snaps.back().t);
    s.plot_time = std::clamp(s.plot_time, t_min, t_max);
    float avail = ImGui::GetContentRegionAvail().x;
    ImGui::SetNextItemWidth(avail - kUnitColWidth);
    ImGui::SliderFloat("##tslider", &s.plot_time, t_min, t_max, "");
    ImGui::SameLine(0, 6);
    ImGui::TextDisabled("s");
    ImGui::Spacing();
    float w = avail - ImGui::CalcTextSize("t =").x - 6.0f - kUnitColWidth;
    ImGui::TextUnformatted("t ="); ImGui::SameLine(0, 6);
    ImGui::SetNextItemWidth(w);
    ExprInputFloat("##texact", &s.plot_time);
    s.plot_time = std::clamp(s.plot_time, t_min, t_max);
    ImGui::SameLine(0, 6); ImGui::TextDisabled("s");
    ImGui::Spacing();
    ImGui::TextDisabled("%d snapshots  [%.4g s … %.4g s]",
                        n_snaps, (double)t_min, (double)t_max);
    ImGui::Spacing();
  }

  if (n_fields > 1) {
    std::vector<const char*> names;
    for (const auto& f : snaps[0].fields) names.push_back(f.name.c_str());
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::Combo("##field", &field_idx, names.data(), n_fields);
    ImGui::Spacing();
  }

  std::vector<double> plot_vals;
  if (transient) {
    plot_vals = InterpolateField(s.result, field_idx, static_cast<double>(s.plot_time));
  } else {
    plot_vals = snaps[0].fields[field_idx].values;
  }

  const std::string& fname = snaps[0].fields[field_idx].name;

  const bool is_electrolyte = is_spme && fname == "c_e";
  const std::vector<double>& x_axis =
      (is_electrolyte && !s.spme_ce_x.empty()) ? s.spme_ce_x : s.cell_centres;
  const char* x_label = is_electrolyte ? "x  [m]"
                      : is_particle    ? "r / R  [-]"
                                       : "x  [m]";
  int n_cells = std::min(static_cast<int>(x_axis.size()),
                         static_cast<int>(plot_vals.size()));

  ImVec2 plot_size = ImGui::GetContentRegionAvail();
  plot_size.y -= 4.0f;

  if (ImPlot::BeginPlot("##results", plot_size)) {
    ImPlot::SetupAxes(x_label, fname.c_str(),
                      ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
    ImPlot::PlotLine(fname.c_str(),
        x_axis.data(), plot_vals.data(), n_cells);
    ImPlot::EndPlot();
  }
}
