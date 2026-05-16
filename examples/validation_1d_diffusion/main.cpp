#include <cmath>
#include <numbers>
#include <print>
#include <string>
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
#include "mphys/transient_solver.hpp"

// ============================================================
// Helpers
// ============================================================

double MaxError(const mphys::Field& numerical, const mphys::Mesh1D& mesh,
                auto analytical_fn) {
  double err = 0.0;
  for (int i = 0; i < mesh.n_cells; ++i) {
    err = std::max(err, std::abs(numerical[i] - analytical_fn(mesh.cell_centres[i])));
  }
  return err;
}

void PrintResult(const std::string& name, double max_error, double tolerance) {
  const bool pass = max_error < tolerance;
  std::println("  {:<45}  L∞ error = {:.2e}   {}",
               name, max_error, pass ? "PASS" : "FAIL");
}

// ============================================================
// Test 1: Steady-state linear diffusion, Dirichlet-Dirichlet
//
// PDE:  D * d²c/dx² = 0   on [0, 1]
// BCs:  c(0) = c_L,  c(1) = c_R
// Analytical:  c(x) = c_L + (c_R - c_L) * x
//
// FVM is exact for linear solutions on any uniform mesh.
// ============================================================

class SteadyLinearModel : public mphys::Model {
 public:
  SteadyLinearModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv,
                    double D, double c_L, double c_R)
      : Model(mesh, sv), D_(D) {
    c_ = AddField("c", 0.5 * (c_L + c_R));  // uniform initial guess
    SetBcs(c_, {mphys::DirichletBc(c_L), mphys::DirichletBc(c_R)});
  }

  void Residual(double /*t*/,
                const std::vector<mphys::Field>& y,
                const std::vector<mphys::Field>& /*ydot*/,
                const std::vector<double>& /*y_alg*/,
                const std::vector<double>& /*ydot_alg*/,
                std::vector<mphys::Field>& rr,
                std::vector<double>& /*rr_alg*/) override {
    // D * d²c/dx² = 0  →  residual = -Laplacian(c, D)
    rr[c_] = -mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_]);
  }

 private:
  int c_ = 0;
  double D_;
};

bool TestSteadyLinear(mphys::SunContext& sunctx, int n_cells) {
  constexpr double kD = 1.0, kCL = 2.0, kCR = 0.5;

  auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, n_cells);
  mphys::StateVector sv(mesh.n_cells);
  SteadyLinearModel model(mesh, sv, kD, kCL, kCR);

  mphys::SolverOptions opts;
  mphys::SteadySolver solver(model, opts, sunctx);
  solver.Solve();

  const double err = MaxError(model.fields()[0], mesh, [](double x) {
    return 2.0 + (0.5 - 2.0) * x;  // c_L + (c_R - c_L)*x
  });
  PrintResult("Steady linear diffusion (n=" + std::to_string(n_cells) + ")", err, 1e-10);
  return err < 1e-10;
}

// ============================================================
// Test 2: Steady-state diffusion with uniform source
//
// PDE:  D * d²c/dx² + S = 0   on [0, 1]
// BCs:  c(0) = 0,  c(1) = 0
// Analytical:  c(x) = S / (2D) * x * (1 - x)
//
// The ghost-cell Dirichlet treatment at boundary cells introduces O(dx²)
// truncation error, giving O(dx²) global accuracy.
// ============================================================

class SteadySourceModel : public mphys::Model {
 public:
  SteadySourceModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv,
                    double D, double S)
      : Model(mesh, sv), D_(D), S_(S) {
    c_ = AddField("c", S / (8.0 * D));  // initial guess ≈ peak of analytical solution
    SetBcs(c_, {mphys::DirichletBc(0.0), mphys::DirichletBc(0.0)});
  }

  void Residual(double /*t*/,
                const std::vector<mphys::Field>& y,
                const std::vector<mphys::Field>& /*ydot*/,
                const std::vector<double>& /*y_alg*/,
                const std::vector<double>& /*ydot_alg*/,
                std::vector<mphys::Field>& rr,
                std::vector<double>& /*rr_alg*/) override {
    // D * d²c/dx² + S = 0  →  residual = -Laplacian(c, D) - S
    mphys::Field lap = mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_]);
    for (int i = 0; i < mesh_.n_cells; ++i) {
      rr[c_][i] = -lap[i] - S_;
    }
  }

 private:
  int c_ = 0;
  double D_;
  double S_;
};

bool TestSteadyWithSource(mphys::SunContext& sunctx, int n_cells) {
  constexpr double kD = 0.5, kS = 2.0;

  auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, n_cells);
  mphys::StateVector sv(mesh.n_cells);
  SteadySourceModel model(mesh, sv, kD, kS);

  mphys::SolverOptions opts;
  mphys::SteadySolver solver(model, opts, sunctx);
  solver.Solve();

  const double err = MaxError(model.fields()[0], mesh, [](double x) {
    return (2.0 / (2.0 * 0.5)) * x * (1.0 - x);  // S/(2D) * x*(1-x)
  });
  // Ghost-cell Dirichlet BCs introduce O(dx²) truncation error at boundary cells.
  const double tolerance = 10.0 / (n_cells * n_cells);
  PrintResult("Steady diffusion + source (n=" + std::to_string(n_cells) + ")", err, tolerance);
  return err < tolerance;
}

// ============================================================
// Test 3: Transient diffusion — decaying sine mode
//
// PDE:  dc/dt = D * d²c/dx²   on [0, 1]
// BCs:  c(0, t) = 0,  c(1, t) = 0
// IC:   c(x, 0) = sin(π*x)
// Analytical:  c(x, t) = sin(π*x) * exp(-D * π² * t)
//
// Spatial error O(dx²); IDA controls temporal accuracy via tolerances.
// ============================================================

class TransientDiffusionModel : public mphys::Model {
 public:
  TransientDiffusionModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv,
                          double D)
      : Model(mesh, sv), D_(D) {
    c_ = AddField("c", 0.0);
    SetBcs(c_, {mphys::DirichletBc(0.0), mphys::DirichletBc(0.0)});
    // Set IC: sin(π*x) at each cell centre
    for (int i = 0; i < mesh.n_cells; ++i) {
      fields_[c_][i] = std::sin(std::numbers::pi * mesh.cell_centres[i]);
    }
  }

  void Residual(double /*t*/,
                const std::vector<mphys::Field>& y,
                const std::vector<mphys::Field>& ydot,
                const std::vector<double>& /*y_alg*/,
                const std::vector<double>& /*ydot_alg*/,
                std::vector<mphys::Field>& rr,
                std::vector<double>& /*rr_alg*/) override {
    // dc/dt - D * d²c/dx² = 0
    rr[c_] = mphys::fvm::Ddt(ydot[c_])
           - mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_]);
  }

 private:
  int c_ = 0;
  double D_;
};

bool TestTransientDiffusion(mphys::SunContext& sunctx, int n_cells, double t_end) {
  constexpr double kD = 1.0;

  auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, n_cells);
  mphys::StateVector sv(mesh.n_cells);
  TransientDiffusionModel model(mesh, sv, kD);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-4;
  opts.maximum_time_step = 0.05;
  opts.tolerance.relative = 1e-8;
  opts.tolerance.absolute = 1e-10;

  mphys::TransientSolver solver(model, opts, sunctx);

  // Only capture the final snapshot
  std::vector<mphys::Field> final_fields;
  solver.Solve(0.0, t_end,
               [&](double /*t*/, const std::vector<mphys::Field>& fields,
                   const std::vector<double>& /*alg*/) {
                 final_fields = fields;
               });

  const double decay = std::exp(-kD * std::numbers::pi * std::numbers::pi * t_end);
  const double err = MaxError(final_fields[0], mesh, [decay](double x) {
    return std::sin(std::numbers::pi * x) * decay;
  });

  // Expected spatial error O(dx²) = O(1/n²)
  const double tolerance = 5.0 / (n_cells * n_cells);
  PrintResult("Transient diffusion t=" + std::to_string(t_end) +
              " (n=" + std::to_string(n_cells) + ")",
              err, tolerance);
  return err < tolerance;
}

// ============================================================
// Mesh refinement study: verify second-order spatial convergence
// ============================================================

void ConvergenceStudy(mphys::SunContext& sunctx) {
  std::println("\n  Mesh refinement study (transient diffusion, t=0.1):");
  std::println("  {:>8}  {:>14}  {:>12}", "n_cells", "L∞ error", "rate");

  constexpr double kD = 1.0, kT = 0.1;
  double prev_err = 0.0;

  for (int n : {25, 50, 100, 200}) {
    auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, n);
    mphys::StateVector sv(mesh.n_cells);
    TransientDiffusionModel model(mesh, sv, kD);

    mphys::SolverOptions opts;
    opts.initial_time_step = 1e-4;
    opts.maximum_time_step = 0.05;
    opts.tolerance.relative = 1e-10;
    opts.tolerance.absolute = 1e-12;

    mphys::TransientSolver solver(model, opts, sunctx);
    std::vector<mphys::Field> final_fields;
    solver.Solve(0.0, kT,
                 [&](double /*t*/, const std::vector<mphys::Field>& f,
                     const std::vector<double>&) { final_fields = f; });

    const double decay = std::exp(-kD * std::numbers::pi * std::numbers::pi * kT);
    const double err = MaxError(final_fields[0], mesh, [decay](double x) {
      return std::sin(std::numbers::pi * x) * decay;
    });

    if (prev_err > 0.0) {
      const double rate = std::log2(prev_err / err);
      std::println("  {:>8}  {:>14.4e}  {:>12.2f}", n, err, rate);
    } else {
      std::println("  {:>8}  {:>14.4e}  {:>12}", n, err, "---");
    }
    prev_err = err;
  }
}

// ============================================================
// main
// ============================================================

int main() {
  std::println("1D Diffusion Validation");
  std::println("{}", std::string(60, '-'));

  mphys::SunContext sunctx;
  int pass = 0, total = 0;

  std::println("\nSteady-state tests:");
  ++total; if (TestSteadyLinear(sunctx, 10))   ++pass;
  ++total; if (TestSteadyLinear(sunctx, 100))  ++pass;
  ++total; if (TestSteadyWithSource(sunctx, 10))  ++pass;
  ++total; if (TestSteadyWithSource(sunctx, 100)) ++pass;

  std::println("\nTransient tests:");
  ++total; if (TestTransientDiffusion(sunctx, 50,  0.1)) ++pass;
  ++total; if (TestTransientDiffusion(sunctx, 100, 0.1)) ++pass;
  ++total; if (TestTransientDiffusion(sunctx, 100, 1.0)) ++pass;

  ConvergenceStudy(sunctx);

  std::println("\n{}", std::string(60, '-'));
  std::println("Result: {}/{} tests passed.", pass, total);

  return pass == total ? 0 : 1;
}
