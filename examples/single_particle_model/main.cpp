#include <cmath>
#include <print>
#include <vector>

#include "mphys/mesh.hpp"
#include "mphys/models/spm.hpp"
#include "mphys/sim_result.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/transient_solver.hpp"

// ============================================================
// Single Particle Model (SPM) — constant-current discharge.
//
// See include/mphys/models/spm.hpp and the PyBaMM SPM notebook:
//   https://docs.pybamm.org/en/stable/source/examples/notebooks/models/SPM.html
//
// Each electrode is one spherical particle solved on a normalised radius
// r ∈ [0,1].  We discharge at constant current and report the terminal-voltage
// curve, then verify lithium couliometry (⟨c_n⟩ falls at exactly -I/(A F ε_n L_n)).
// ============================================================

using mphys::models::SpmModel;
using mphys::models::SpmParameters;

static double SphericalAverage(const mphys::Field& c, const mphys::Mesh1D& mesh) {
  double num = 0.0, den = 0.0;
  for (int i = 0; i < mesh.n_cells; ++i) {
    const double rl = mesh.face_positions[i];
    const double rr = mesh.face_positions[i + 1];
    const double vol = (rr * rr * rr - rl * rl * rl) / 3.0;
    num += vol * c[i];
    den += vol;
  }
  return num / den;
}

int main() {
  SpmParameters p;          // LG M50 (Chen 2020) defaults
  p.I = 5.0;                // ~1C discharge [A]

  constexpr int    kNCells = 30;
  constexpr double kTEnd   = 3000.0;  // [s]

  auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, kNCells,
                                       mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(mesh.n_cells);
  SpmModel model(mesh, sv, p);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-2;
  opts.maximum_time_step = 20.0;
  opts.tolerance.relative = 1e-7;
  opts.tolerance.absolute = 1e-7;

  mphys::SunContext sunctx;
  mphys::SimResult  result;

  mphys::TransientSolver solver(model, opts, sunctx);
  const std::string warn = solver.Solve(
      0.0, kTEnd,
      [&](double t, const std::vector<mphys::Field>& f,
          const std::vector<double>& a) { result.Record(t, f, a); });

  std::println("Single Particle Model — {:.1f} A discharge", p.I);
  std::println("{}", std::string(60, '-'));
  if (!warn.empty()) std::println("solver: {}", warn);
  std::println("{} snapshots over {:.0f} s\n", result.snapshots.size(), kTEnd);

  std::println("  {:>10}  {:>12}  {:>14}  {:>14}",
               "t [s]", "V [V]", "c_n,surf", "c_p,surf");
  const int n = static_cast<int>(result.snapshots.size());
  for (int s = 0; s < n; s += std::max(1, n / 12)) {
    const auto& snap = result.snapshots[s];
    const auto& cn = snap.fields[model.c_n_index()];
    const auto& cp = snap.fields[model.c_p_index()];
    std::println("  {:>10.1f}  {:>12.4f}  {:>14.1f}  {:>14.1f}",
                 snap.t, snap.algebraics[model.voltage_index()],
                 cn[kNCells - 1], cp[kNCells - 1]);
  }

  // Couliometric check: ⟨c_n⟩(t) - c_n0 should equal -I/(A F ε_n L_n) · t.
  const auto& final_cn = result.snapshots.back().fields[model.c_n_index()];
  const double t_final = result.snapshots.back().t;
  const double avg = SphericalAverage(final_cn, mesh);
  const double rate = -p.I / (p.A * SpmParameters::F * p.eps_n * p.L_n);
  const double expected = p.cn0 + rate * t_final;
  const double err = std::abs(avg - expected);

  std::println("\nCouliometry check at t = {:.0f} s:", t_final);
  std::println("  <c_n>      = {:.4f} mol/m^3", avg);
  std::println("  expected   = {:.4f} mol/m^3", expected);
  std::println("  |error|    = {:.3e}", err);

  const bool pass = err < 1e-2;
  std::println("\nResult: {}", pass ? "PASS" : "FAIL");
  return pass ? 0 : 1;
}
