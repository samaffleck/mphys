#include <cmath>
#include <print>
#include <vector>

#include "mphys/mesh.hpp"
#include "mphys/models/spme.hpp"
#include "mphys/sim_result.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/transient_solver.hpp"

// ============================================================
// Single Particle Model with Electrolyte (SPMe) — constant-current discharge.
//
// See include/mphys/models/spme.hpp and the PyBaMM SPMe notebook:
//   https://docs.pybamm.org/en/stable/source/examples/notebooks/models/SPMe.html
//
// Two spherical particles (normalised r ∈ [0,1]) plus a macroscopic electrolyte
// concentration on the negative|separator|positive sandwich.  We discharge at
// constant current, report the terminal-voltage curve, and verify both lithium
// couliometry in the particle and porosity-weighted electrolyte conservation.
// ============================================================

using mphys::models::MakeSpmeElectrolyteMesh;
using mphys::models::SpmeModel;
using mphys::models::SpmeParameters;
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
  SpmeParameters p;          // LG M50 (Chen 2020) defaults
  p.core.I = 5.0;            // ~1C discharge [A]

  constexpr int    kNn = 25, kNs = 15, kNp = 25;
  constexpr double kTEnd = 3000.0;  // [s]

  auto em = MakeSpmeElectrolyteMesh(p, kNn, kNs, kNp);
  auto particle_mesh = mphys::MakeUniformMesh1D(
      0.0, 1.0, em.mesh.n_cells, mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(em.mesh.n_cells);
  SpmeModel model(particle_mesh, em, sv, p);

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

  std::println("Single Particle Model w/ Electrolyte — {:.1f} A discharge", p.core.I);
  std::println("{}", std::string(72, '-'));
  if (!warn.empty()) std::println("solver: {}", warn);
  std::println("{} snapshots over {:.0f} s\n", result.snapshots.size(), kTEnd);

  std::println("  {:>10}  {:>12}  {:>14}  {:>14}  {:>12}",
               "t [s]", "V [V]", "c_n,surf", "c_p,surf", "c_e (sep)");
  const int n = static_cast<int>(result.snapshots.size());
  const int sep_cell = kNn + kNs / 2;   // a separator cell for monitoring
  for (int s = 0; s < n; s += std::max(1, n / 12)) {
    const auto& snap = result.snapshots[s];
    const auto& cn = snap.fields[model.c_n_index()];
    const auto& cp = snap.fields[model.c_p_index()];
    const auto& ce = snap.fields[model.c_e_index()];
    std::println("  {:>10.1f}  {:>12.4f}  {:>14.1f}  {:>14.1f}  {:>12.2f}",
                 snap.t, snap.algebraics[model.voltage_index()],
                 cn[em.mesh.n_cells - 1], cp[em.mesh.n_cells - 1], ce[sep_cell]);
  }

  // Couliometry: ⟨c_n⟩(t) - c_n0 should equal -I/(A F ε_n L_n) · t.
  const auto& final_cn = result.snapshots.back().fields[model.c_n_index()];
  const double t_final = result.snapshots.back().t;
  const double avg = SphericalAverage(final_cn, particle_mesh);
  const double rate = -p.core.I /
                      (p.core.A * SpmParameters::F * p.core.eps_n * p.core.L_n);
  const double expected = p.core.cn0 + rate * t_final;
  const double err_coul = std::abs(avg - expected);

  // Electrolyte conservation: ⟨ε c_e⟩ must be unchanged.
  const auto& final_ce = result.snapshots.back().fields[model.c_e_index()];
  auto eps_mass = [&](const auto& ce) {
    double m = 0.0;
    for (int i = 0; i < em.mesh.n_cells; ++i) {
      const double eps = em.region[i] == 0 ? p.eps_e_n
                       : em.region[i] == 1 ? p.eps_e_s : p.eps_e_p;
      m += eps * ce[i] * em.mesh.dx[i];
    }
    return m;
  };
  const mphys::Field ce0("c_e", em.mesh.n_cells, p.ce0);
  const double mass0 = eps_mass(ce0);
  const double mass1 = eps_mass(final_ce);
  const double err_mass = std::abs(mass1 - mass0) / mass0;

  std::println("\nCouliometry check at t = {:.0f} s:", t_final);
  std::println("  <c_n>     = {:.4f} mol/m^3", avg);
  std::println("  expected  = {:.4f} mol/m^3", expected);
  std::println("  |error|   = {:.3e}", err_coul);
  std::println("\nElectrolyte conservation:");
  std::println("  <eps c_e> drift (relative) = {:.3e}", err_mass);

  const bool pass = err_coul < 1e-2 && err_mass < 1e-6;
  std::println("\nResult: {}", pass ? "PASS" : "FAIL");
  return pass ? 0 : 1;
}
