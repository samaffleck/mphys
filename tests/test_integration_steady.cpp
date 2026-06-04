// Integration tests: full SteadySolver (KINSOL) runs against models with known
// analytical solutions. These exercise StateVector, Model, the FVM operators
// and the SUNDIALS wrapper end to end.
#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/fvm_operators.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/steady_solver.hpp"
#include "mphys/sun_context.hpp"

namespace {

using mphys::DirichletBc;
using mphys::Field;
using mphys::MakeUniformMesh1D;
using mphys::Mesh1D;
using mphys::Model;
using mphys::StateVector;

double MaxError(const Field& numerical, const Mesh1D& mesh, auto analytical) {
  double err = 0.0;
  for (int i = 0; i < mesh.n_cells; ++i)
    err = std::max(err, std::abs(numerical[i] - analytical(mesh.cell_centres[i])));
  return err;
}

// D c'' = 0,  c(0)=c_L, c(1)=c_R  ->  c(x) = c_L + (c_R - c_L) x
class SteadyLinearModel : public Model {
 public:
  SteadyLinearModel(const Mesh1D& mesh, StateVector& sv, double D, double cL,
                    double cR)
      : Model(mesh, sv), D_(D) {
    c_ = AddField("c", 0.5 * (cL + cR));
    SetBcs(c_, {DirichletBc(cL), DirichletBc(cR)});
  }
  void Residual(double, const std::vector<Field>& y, const std::vector<Field>&,
                const std::vector<double>&, const std::vector<double>&,
                std::vector<Field>& rr, std::vector<double>&) override {
    rr[c_] = -mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_]);
  }

 private:
  int c_ = 0;
  double D_;
};

// D c'' + S = 0, c(0)=c(1)=0  ->  c(x) = S/(2D) x (1-x)
class SteadySourceModel : public Model {
 public:
  SteadySourceModel(const Mesh1D& mesh, StateVector& sv, double D, double S)
      : Model(mesh, sv), D_(D), S_(S) {
    c_ = AddField("c", S / (8.0 * D));
    SetBcs(c_, {DirichletBc(0.0), DirichletBc(0.0)});
  }
  void Residual(double, const std::vector<Field>& y, const std::vector<Field>&,
                const std::vector<double>&, const std::vector<double>&,
                std::vector<Field>& rr, std::vector<double>&) override {
    Field lap = mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_]);
    for (int i = 0; i < mesh_.n_cells; ++i) rr[c_][i] = -lap[i] - S_;
  }

 private:
  int c_ = 0;
  double D_, S_;
};

// Spherical: (1/r²)(r² D c')' = 0, c(r_i)=c_i, c(r_o)=c_o  ->  c(r)=A+B/r
class SphericalDiffusionModel : public Model {
 public:
  SphericalDiffusionModel(const Mesh1D& mesh, StateVector& sv, double D,
                          double c_inner, double c_outer)
      : Model(mesh, sv), D_(D) {
    c_ = AddField("c", 0.5 * (c_inner + c_outer));
    SetBcs(c_, {DirichletBc(c_inner), DirichletBc(c_outer)});
  }
  void Residual(double, const std::vector<Field>& y, const std::vector<Field>&,
                const std::vector<double>&, const std::vector<double>&,
                std::vector<Field>& rr, std::vector<double>&) override {
    rr[c_] = -mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_]);
  }

 private:
  int c_ = 0;
  double D_;
};

TEST(SteadyIntegration, LinearDiffusionIsExact) {
  const double D = 1.0, cL = 2.0, cR = 0.5;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 20);
  StateVector sv(mesh.n_cells);
  SteadyLinearModel model(mesh, sv, D, cL, cR);

  mphys::SunContext ctx;
  mphys::SolverOptions opts;
  mphys::SteadySolver solver(model, opts, ctx);
  solver.Solve();

  // FVM is exact for linear solutions regardless of mesh resolution.
  const double err =
      MaxError(model.fields()[0], mesh, [&](double x) { return cL + (cR - cL) * x; });
  EXPECT_LT(err, 1e-10);
}

TEST(SteadyIntegration, DiffusionWithSourceConvergesSecondOrder) {
  const double D = 0.5, S = 2.0;
  auto run = [&](int n) {
    Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, n);
    StateVector sv(mesh.n_cells);
    SteadySourceModel model(mesh, sv, D, S);
    mphys::SunContext ctx;
    mphys::SolverOptions opts;
    mphys::SteadySolver solver(model, opts, ctx);
    solver.Solve();
    return MaxError(model.fields()[0], mesh,
                    [&](double x) { return S / (2.0 * D) * x * (1.0 - x); });
  };

  const double err_coarse = run(20);
  const double err_fine = run(40);
  EXPECT_LT(err_fine, 10.0 / (40 * 40));
  // Halving dx should cut the error by ~4x (second order). Allow slack.
  EXPECT_GT(err_coarse / err_fine, 3.0);
}

TEST(SteadyIntegration, LineSearchStrategyAlsoConverges) {
  const double D = 1.0, cL = 2.0, cR = 0.5;
  Mesh1D mesh = MakeUniformMesh1D(0.0, 1.0, 20);
  StateVector sv(mesh.n_cells);
  SteadyLinearModel model(mesh, sv, D, cL, cR);

  mphys::SunContext ctx;
  mphys::SolverOptions opts;
  mphys::SteadySolver solver(model, opts, ctx);
  solver.Solve(mphys::NewtonStrategy::kLineSearch);

  const double err =
      MaxError(model.fields()[0], mesh, [&](double x) { return cL + (cR - cL) * x; });
  EXPECT_LT(err, 1e-8);
}

TEST(SteadyIntegration, SphericalDiffusionMatchesAnalytical) {
  const double D = 1.0, r_i = 1.0, r_o = 2.0, c_i = 2.0, c_o = 0.5;
  const double B = (c_i - c_o) / (1.0 / r_i - 1.0 / r_o);
  const double A = c_i - B / r_i;

  const int n = 100;
  Mesh1D mesh = MakeUniformMesh1D(r_i, r_o, n, mphys::CoordSystem::kSpherical);
  StateVector sv(mesh.n_cells);
  SphericalDiffusionModel model(mesh, sv, D, c_i, c_o);

  mphys::SunContext ctx;
  mphys::SolverOptions opts;
  mphys::SteadySolver solver(model, opts, ctx);
  solver.Solve();

  const double err =
      MaxError(model.fields()[0], mesh, [&](double r) { return A + B / r; });
  const double dr = (r_o - r_i) / n;
  EXPECT_LT(err, 10.0 * dr * dr);
}

}  // namespace
