#include <cmath>
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

// ============================================================
// Steady-state spherical diffusion
//
// PDE:  (1/r²) d/dr(r² D dc/dr) = 0   on [r_inner, r_outer]
// BCs:  c(r_inner) = c_inner  (Dirichlet)
//       c(r_outer) = c_outer  (Dirichlet)
//
// Analytical solution:  c(r) = A + B/r
//   B = (c_inner - c_outer) / (1/r_inner - 1/r_outer)
//   A = c_inner - B/r_inner
//
// The FVM stencil reduces to the correct spherical form once
// CoordSystem::kSpherical is passed to MakeUniformMesh1D.
// ============================================================

class SphericalDiffusionModel : public mphys::Model {
 public:
  SphericalDiffusionModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv,
                          double D, double c_inner, double c_outer)
      : Model(mesh, sv), D_(D) {
    c_ = AddField("c", 0.5 * (c_inner + c_outer));
    SetBcs(c_, {mphys::DirichletBc(c_inner), mphys::DirichletBc(c_outer)});
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

// Analytical solution coefficients
static void AnalyticalCoeffs(double r_inner, double r_outer,
                              double c_inner, double c_outer,
                              double& A, double& B) {
  B = (c_inner - c_outer) / (1.0 / r_inner - 1.0 / r_outer);
  A = c_inner - B / r_inner;
}

static double Analytical(double r, double A, double B) {
  return A + B / r;
}

static double MaxError(const mphys::Field& numerical, const mphys::Mesh1D& mesh,
                       double A, double B) {
  double err = 0.0;
  for (int i = 0; i < mesh.n_cells; ++i) {
    err = std::max(err, std::abs(numerical[i] - Analytical(mesh.cell_centres[i], A, B)));
  }
  return err;
}

static bool RunTest(mphys::SunContext& sunctx, int n_cells,
                    double r_inner, double r_outer,
                    double c_inner, double c_outer, double D) {
  auto mesh = mphys::MakeUniformMesh1D(r_inner, r_outer, n_cells,
                                        mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(mesh.n_cells);
  SphericalDiffusionModel model(mesh, sv, D, c_inner, c_outer);

  mphys::SolverOptions opts;
  mphys::SteadySolver solver(model, opts, sunctx);
  solver.Solve();

  double A, B;
  AnalyticalCoeffs(r_inner, r_outer, c_inner, c_outer, A, B);
  const double err = MaxError(model.fields()[0], mesh, A, B);

  // Second-order spatial convergence: error ≈ C/n²
  const double dr = (r_outer - r_inner) / n_cells;
  const double tol = 10.0 * dr * dr;
  const bool pass = err < tol;
  std::println("  n={:4d}  L∞ err = {:.3e}  tol = {:.3e}  {}",
               n_cells, err, tol, pass ? "PASS" : "FAIL");
  return pass;
}

int main() {
  std::println("Spherical Diffusion Validation");
  std::println("{}", std::string(60, '-'));

  mphys::SunContext sunctx;

  constexpr double kD      = 1.0;
  constexpr double kRInner = 1.0;
  constexpr double kROuter = 2.0;
  constexpr double kCInner = 2.0;
  constexpr double kCOuter = 0.5;

  double A, B;
  AnalyticalCoeffs(kRInner, kROuter, kCInner, kCOuter, A, B);
  std::println("Analytical: c(r) = {:.4f} + {:.4f}/r", A, B);
  std::println("");

  int pass = 0, total = 0;
  for (int n : {10, 20, 50, 100, 200}) {
    ++total;
    if (RunTest(sunctx, n, kRInner, kROuter, kCInner, kCOuter, kD)) ++pass;
  }

  // Convergence rate study
  std::println("\n  Mesh refinement study:");
  std::println("  {:>6}  {:>14}  {:>10}", "n_cells", "L∞ error", "rate");
  double prev_err = 0.0;
  for (int n : {25, 50, 100, 200}) {
    auto mesh = mphys::MakeUniformMesh1D(kRInner, kROuter, n,
                                          mphys::CoordSystem::kSpherical);
    mphys::StateVector sv(mesh.n_cells);
    SphericalDiffusionModel model(mesh, sv, kD, kCInner, kCOuter);
    mphys::SolverOptions opts;
    mphys::SteadySolver solver(model, opts, sunctx);
    solver.Solve();

    double A2, B2;
    AnalyticalCoeffs(kRInner, kROuter, kCInner, kCOuter, A2, B2);
    const double err = MaxError(model.fields()[0], mesh, A2, B2);

    if (prev_err > 0.0) {
      std::println("  {:>6}  {:>14.4e}  {:>10.2f}", n, err, std::log2(prev_err / err));
    } else {
      std::println("  {:>6}  {:>14.4e}  {:>10}", n, err, "---");
    }
    prev_err = err;
  }

  std::println("\n{}", std::string(60, '-'));
  std::println("Result: {}/{} tests passed.", pass, total);
  return pass == total ? 0 : 1;
}
