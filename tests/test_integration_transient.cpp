// Integration tests: full TransientSolver (IDA) runs against models with known
// analytical solutions, plus exercises the output-callback mechanism.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/fvm_operators.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/sim_result.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/transient_solver.hpp"

namespace {

using mphys::DirichletBc;
using mphys::Field;
using mphys::MakeUniformMesh1D;
using mphys::Mesh1D;
using mphys::Model;
using mphys::NeumannBc;
using mphys::StateVector;

double MaxError(const Field& numerical, const Mesh1D& mesh, auto analytical) {
  double err = 0.0;
  for (int i = 0; i < mesh.n_cells; ++i)
    err = std::max(err, std::abs(numerical[i] - analytical(mesh.cell_centres[i])));
  return err;
}

// dc/dt = D c'', c(0)=c(1)=0, c(x,0)=sin(pi x)
// Analytical: c(x,t) = sin(pi x) exp(-D pi^2 t)
class TransientDiffusionModel : public Model {
 public:
  TransientDiffusionModel(const Mesh1D& mesh, StateVector& sv, double D)
      : Model(mesh, sv), D_(D) {
    c_ = AddField("c", 0.0);
    SetBcs(c_, {DirichletBc(0.0), DirichletBc(0.0)});
    for (int i = 0; i < mesh.n_cells; ++i)
      fields_[c_][i] = std::sin(std::numbers::pi * mesh.cell_centres[i]);
  }
  void Residual(double, const std::vector<Field>& y, const std::vector<Field>& ydot,
                const std::vector<double>&, const std::vector<double>&,
                std::vector<Field>& rr, std::vector<double>&) override {
    rr[c_] = mphys::fvm::Ddt(ydot[c_]) -
             mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_]);
  }

 private:
  int c_ = 0;
  double D_;
};

// Convection-diffusion-reaction in a packed bed; Dirichlet inlet, Neumann outlet.
class ReactorModel : public Model {
 public:
  ReactorModel(const Mesh1D& mesh, StateVector& sv, double D, double u, double k)
      : Model(mesh, sv), D_(D), u_(u), k_(k) {
    c_ = AddField("c", 0.0);
    SetBcs(c_, {DirichletBc(1.0), NeumannBc(0.0)});
  }
  void Residual(double, const std::vector<Field>& y, const std::vector<Field>& ydot,
                const std::vector<double>&, const std::vector<double>&,
                std::vector<Field>& rr, std::vector<double>&) override {
    rr[c_] = mphys::fvm::Ddt(ydot[c_]) +
             mphys::fvm::Convection(y[c_], u_, mesh_, bcs_[c_]) -
             mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_]) + y[c_] * k_;
  }

 private:
  int c_ = 0;
  double D_, u_, k_;
};

mphys::SolverOptions AccurateOptions() {
  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-4;
  opts.maximum_time_step = 0.05;
  opts.tolerance.relative = 1e-8;
  opts.tolerance.absolute = 1e-10;
  return opts;
}

TEST(TransientIntegration, DecayingSineModeMatchesAnalytical) {
  const double D = 1.0, t_end = 0.1;
  const int n = 100;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, n);
  StateVector sv(mesh.n_cells);
  TransientDiffusionModel model(mesh, sv, D);

  mphys::SunContext ctx;
  mphys::TransientSolver solver(model, AccurateOptions(), ctx);

  std::vector<Field> final_fields;
  const std::string msg = solver.Solve(
      0.0, t_end,
      [&](double, const std::vector<Field>& f, const std::vector<double>&) {
        final_fields = f;
      });

  EXPECT_TRUE(msg.empty()) << msg;
  const double decay = std::exp(-D * std::numbers::pi * std::numbers::pi * t_end);
  const double err = MaxError(final_fields[0], mesh, [&](double x) {
    return std::sin(std::numbers::pi * x) * decay;
  });
  EXPECT_LT(err, 5.0 / (n * n));
}

TEST(TransientIntegration, SolutionDecaysMonotonically) {
  const double D = 1.0;
  const int n = 50;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, n);
  StateVector sv(mesh.n_cells);
  TransientDiffusionModel model(mesh, sv, D);

  mphys::SunContext ctx;
  mphys::TransientSolver solver(model, AccurateOptions(), ctx);

  // Peak amplitude (centre cell) must decay over time.
  const int mid = n / 2;
  double prev_peak = 1e9;
  bool first = true;
  std::string msg = solver.Solve(
      0.0, 0.2,
      [&](double, const std::vector<Field>& f, const std::vector<double>&) {
        const double peak = f[0][mid];
        if (!first) EXPECT_LE(peak, prev_peak + 1e-9);
        prev_peak = peak;
        first = false;
      });
  EXPECT_TRUE(msg.empty()) << msg;
}

TEST(TransientIntegration, CallbackRecordsInitialStateAndAdvancingTime) {
  const int n = 40;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, n);
  StateVector sv(mesh.n_cells);
  ReactorModel model(mesh, sv, 1e-4, 0.01, 0.1);

  mphys::SunContext ctx;
  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-3;
  opts.maximum_time_step = 1.0;
  opts.tolerance.relative = 1e-6;
  opts.tolerance.absolute = 1e-8;
  mphys::TransientSolver solver(model, opts, ctx);

  mphys::SimResult result;
  const double t_end = 10.0;
  const std::string msg = solver.Solve(
      0.0, t_end,
      [&](double t, const std::vector<Field>& f, const std::vector<double>& a) {
        result.Record(t, f, a);
      });
  EXPECT_TRUE(msg.empty()) << msg;

  ASSERT_GE(result.snapshots.size(), 2u);
  // First snapshot is the initial condition at t0.
  EXPECT_DOUBLE_EQ(result.snapshots.front().t, 0.0);
  // Times are non-decreasing and reach the end.
  for (size_t i = 1; i < result.snapshots.size(); ++i)
    EXPECT_GE(result.snapshots[i].t, result.snapshots[i - 1].t);
  EXPECT_NEAR(result.snapshots.back().t, t_end, 1e-6);
}

TEST(TransientIntegration, DirichletInletBoundaryIsRespected) {
  const int n = 50;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, n);
  StateVector sv(mesh.n_cells);
  ReactorModel model(mesh, sv, 1e-3, 0.05, 0.0);

  mphys::SunContext ctx;
  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-3;
  opts.maximum_time_step = 1.0;
  opts.tolerance.relative = 1e-6;
  opts.tolerance.absolute = 1e-8;
  mphys::TransientSolver solver(model, opts, ctx);

  std::vector<Field> final_fields;
  solver.Solve(0.0, 100.0,
               [&](double, const std::vector<Field>& f, const std::vector<double>&) {
                 final_fields = f;
               });

  // With inlet c(0)=1 held and Neumann outlet, the near-inlet cell sits between
  // the inlet value and the interior, and all values stay within [0, 1].
  for (int i = 0; i < n; ++i) {
    EXPECT_GE(final_fields[0][i], -1e-6);
    EXPECT_LE(final_fields[0][i], 1.0 + 1e-6);
  }
  // At steady state with no reaction the first cell is close to the inlet value.
  EXPECT_GT(final_fields[0][0], 0.5);
}

}  // namespace
