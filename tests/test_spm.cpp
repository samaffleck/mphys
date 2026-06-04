#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "mphys/mesh.hpp"
#include "mphys/models/spm.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/transient_solver.hpp"

namespace {

using mphys::models::SpmModel;
using mphys::models::SpmParameters;

// Volume-weighted average of a cell-centred field over the spherical mesh.
double SphericalAverage(const mphys::Field& c, const mphys::Mesh1D& mesh) {
  double num = 0.0, den = 0.0;
  for (int i = 0; i < mesh.n_cells; ++i) {
    const double rl = mesh.face_positions[i];
    const double rr = mesh.face_positions[i + 1];
    const double vol = (rr * rr * rr - rl * rl * rl) / 3.0;  // 4π/3 cancels
    num += vol * c[i];
    den += vol;
  }
  return num / den;
}

// Build SPM parameters tuned so the negative particle reduces to the canonical
// unit-sphere diffusion problem (R=D=ε=L=A=1, i_app = F):
//   D_eff,n = 1,  surface gradient g_n = -1/3,  d⟨c_n⟩/dt = -1.
// This makes the expected values exact and easy to verify by hand.
SpmParameters MakeTestParams() {
  SpmParameters p;
  p.R_n = 1.0;  p.D_n = 1.0;  p.eps_n = 1.0;  p.L_n = 1.0;
  p.A = 1.0;    p.I = SpmParameters::F;          // i_app = I/A = F
  p.cn0 = 100.0; p.cn_max = 1000.0;
  // Positive particle: any physically valid configuration (not under test).
  p.R_p = 1.0;  p.D_p = 1.0;  p.eps_p = 1.0;  p.L_p = 1.0;
  p.cp0 = 300.0; p.cp_max = 1000.0;
  p.c_e = 1000.0;
  return p;
}

}  // namespace

// Couliometry: the volume-averaged negative concentration must fall at exactly
// the rate dictated by the applied current, d⟨c_n⟩/dt = -I/(A F ε_n L_n).
// With the test parameters this rate is -1, so ⟨c_n⟩(t_end=1) = c_n0 - 1 = 99.
TEST(Spm, LithiumConservation) {
  const SpmParameters p = MakeTestParams();
  const double t_end = 1.0;

  auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, 80, mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(mesh.n_cells);
  SpmModel model(mesh, sv, p);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-4;
  opts.maximum_time_step = 0.05;
  opts.tolerance.relative = 1e-10;
  opts.tolerance.absolute = 1e-10;

  mphys::SunContext sunctx;
  mphys::TransientSolver solver(model, opts, sunctx);

  std::vector<mphys::Field> final_fields;
  solver.Solve(0.0, t_end,
               [&](double, const std::vector<mphys::Field>& f,
                   const std::vector<double>&) { final_fields = f; });

  const double rate = -p.I / (p.A * SpmParameters::F * p.eps_n * p.L_n);  // = -1
  const double expected_avg = p.cn0 + rate * t_end;                      // = 99
  const double avg = SphericalAverage(final_fields[model.c_n_index()], mesh);

  EXPECT_NEAR(avg, expected_avg, 1e-4);
}

// Constant-current discharge drives a constant surface flux, so at large times
// the negative particle relaxes onto the Carslaw & Jaeger quasi-steady profile
// (Crank, "The Mathematics of Diffusion", §6.3):
//   c(r,t) = c0 + 3 D_eff g t + g ( r²/2 - 3/10 )       (a = 1)
// where g is the surface gradient (Neumann BC).  With the test parameters
// D_eff = 1, g = -1/3, t = 1  ⇒  c(r,1) = 99.1 - r²/6.
TEST(Spm, CarslawJaegerProfile) {
  const SpmParameters p = MakeTestParams();
  const double t_end = 1.0;
  const double D_eff = p.Dn_eff();                                  // = 1
  const double g = -p.R_n * p.R_n * p.i_app() /
                   (3.0 * p.D_n * SpmParameters::F * p.eps_n * p.L_n);  // = -1/3

  auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, 120, mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(mesh.n_cells);
  SpmModel model(mesh, sv, p);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-4;
  opts.maximum_time_step = 0.02;
  opts.tolerance.relative = 1e-10;
  opts.tolerance.absolute = 1e-10;

  mphys::SunContext sunctx;
  mphys::TransientSolver solver(model, opts, sunctx);

  std::vector<mphys::Field> final_fields;
  solver.Solve(0.0, t_end,
               [&](double, const std::vector<mphys::Field>& f,
                   const std::vector<double>&) { final_fields = f; });

  const auto& cn = final_fields[model.c_n_index()];
  double max_err = 0.0;
  for (int i = 0; i < mesh.n_cells; ++i) {
    const double r = mesh.cell_centres[i];
    const double analytic = p.cn0 + 3.0 * D_eff * g * t_end
                          + g * (0.5 * r * r - 0.3);
    max_err = std::max(max_err, std::abs(cn[i] - analytic));
  }

  // Second-order spatial error plus a negligible decayed transient (~e^-20).
  EXPECT_LT(max_err, 2e-3);
}

// The reported voltage must equal the closed-form SPM expression evaluated at
// the surface concentrations, and sit in a physically sensible Li-ion range.
TEST(Spm, VoltageOutput) {
  SpmParameters p;            // realistic LG M50 defaults
  p.I = 5.0;
  const double t_end = 50.0;

  auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, 30, mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(mesh.n_cells);
  SpmModel model(mesh, sv, p);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-3;
  opts.maximum_time_step = 5.0;
  opts.tolerance.relative = 1e-8;
  opts.tolerance.absolute = 1e-8;

  mphys::SunContext sunctx;
  mphys::TransientSolver solver(model, opts, sunctx);

  std::vector<mphys::Field> ff;
  std::vector<double> fa;
  solver.Solve(0.0, t_end,
               [&](double, const std::vector<mphys::Field>& f,
                   const std::vector<double>& a) { ff = f; fa = a; });

  // Surface concentrations via the model's own extrapolation (gradient signs
  // match the Neumann BCs the model installs internally).
  const double gn = -p.R_n * p.R_n * p.i_app() /
                    (3.0 * p.D_n * SpmParameters::F * p.eps_n * p.L_n);
  const double gp = +p.R_p * p.R_p * p.i_app() /
                    (3.0 * p.D_p * SpmParameters::F * p.eps_p * p.L_p);
  const double cn_surf = model.SurfaceValue(ff[model.c_n_index()], gn);
  const double cp_surf = model.SurfaceValue(ff[model.c_p_index()], gp);

  const double v_expected = mphys::models::SpmVoltage(p, cn_surf, cp_surf);
  const double v_solver = fa[model.voltage_index()];

  EXPECT_NEAR(v_solver, v_expected, 1e-6);
  EXPECT_GT(v_solver, 2.0);
  EXPECT_LT(v_solver, 4.5);
}
