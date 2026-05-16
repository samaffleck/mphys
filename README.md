# mphys

A C++23 finite volume library for solving 1D PDEs. Physics models are defined by subclassing `Model` and implementing a residual function; the library handles mesh generation, state management, and time integration via SUNDIALS.

## Features

- Cell-centred FVM on 1D structured meshes — Cartesian, cylindrical, and spherical coordinates
- Transient integration (SUNDIALS IDA) and steady-state solving (SUNDIALS KINSOL)
- FVM operators: `Laplacian`, `Convection`, `Grad`, `Div`, `Ddt`
- COMSOL-inspired desktop GUI with live results plotting (ImGui + ImPlot)

## Requirements

- CMake ≥ 3.25, Ninja
- Clang with C++23 support (`/opt/homebrew/opt/llvm/bin/clang++`)
- GLFW (`brew install glfw`) — GUI only

All other dependencies (SUNDIALS, Google Test, ImGui, ImPlot) are included as git submodules.

```bash
git clone --recurse-submodules <url>
```

## Build

```bash
cmake --preset mac-debug    # configure
cmake --build --preset mac-debug

cmake --preset mac-release
cmake --build --preset mac-release
```

Targets of interest:

| Target | Path |
|--------|------|
| `mphys_gui` | `build/mac-debug/gui/mphys_gui` |
| `mphys_tests` | `build/mac-debug/tests/mphys_tests` |
| `example_convection_diffusion` | `build/mac-debug/examples/...` |
| `example_spherical_diffusion` | `build/mac-debug/examples/...` |

To build the GUI only (faster iteration):

```bash
cmake --build --preset mac-debug-gui
```

## Defining a physics model

Subclass `Model`, declare fields in the constructor, and implement `Residual`. The same implementation is used for both transient and steady-state solves — `ydot` is empty for steady-state.

```cpp
class ReactorModel : public mphys::Model {
 public:
  ReactorModel(const mphys::Mesh1D& mesh, mphys::StateVector& sv, double D, double u)
      : Model(mesh, sv), D_(D), u_(u) {
    c_ = AddField("c", 0.0);
    SetBcs(c_, {mphys::DirichletBc(1.0), mphys::NeumannBc(0.0)});
  }

  void Residual(double, const std::vector<mphys::Field>& y,
                const std::vector<mphys::Field>& ydot,
                const std::vector<double>&, const std::vector<double>&,
                std::vector<mphys::Field>& rr, std::vector<double>&) override {
    rr[c_] = mphys::fvm::Ddt(ydot[c_])
           + mphys::fvm::Convection(y[c_], u_, mesh_, bcs_[c_])
           - mphys::fvm::Laplacian(y[c_], D_, mesh_, bcs_[c_]);
  }

 private:
  int c_ = 0;
  double D_, u_;
};
```

Running a transient solve:

```cpp
auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, 100);
mphys::StateVector sv(mesh.n_cells);
ReactorModel model(mesh, sv, 1e-4, 0.01);

mphys::SolverOptions opts;
mphys::SunContext sunctx;
mphys::SimResult result;

mphys::TransientSolver solver(model, opts, sunctx);
solver.Solve(0.0, 50.0, [&](double t, const auto& fields, const auto& alg) {
    result.Record(t, fields, alg);
});
```

Spherical coordinates require only a different mesh factory call — the `Residual` is identical:

```cpp
auto mesh = mphys::MakeUniformMesh1D(1.0, 2.0, 100, mphys::CoordSystem::kSpherical);
```

## License

For academic and non-commercial use only. © 2026 Sam Affleck and contributors.
