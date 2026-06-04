# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

CMake presets are defined in `CMakePresets.json` (uses Ninja + clang from Homebrew LLVM). Both presets enable `BUILD_GUI=ON`.

```bash
# Configure
cmake --preset mac-debug
cmake --preset mac-release

# Build everything
cmake --build --preset mac-debug
cmake --build --preset mac-release

# Build GUI only (fastest iteration on GUI changes)
cmake --build --preset mac-debug-gui
cmake --build --preset mac-release-gui
```

VS Code launch configs ("mphys GUI (Debug/Release)", etc.) and matching build tasks are in `.vscode/`.

## Tests

```bash
# Build and run all tests
cmake --build build/mac-debug --target mphys_tests -j8
cd build/mac-debug && ctest -L mphys --output-on-failure   # -L mphys excludes SUNDIALS' own tests

# Or run the gtest binary directly (faster, full gtest output)
./build/mac-debug/tests/mphys_tests

# Run the validation suite (analytical convergence checks — more informative than gtest)
./build/mac-debug/examples/example_validation_1d_diffusion
```

Unit tests live in `tests/test_<unit>.cpp` (one file per translation unit:
mesh, field, boundary_condition, state_vector, fvm_operators, model,
solver_options, sundials_types). Integration tests that drive the full solvers
against analytical solutions are in `tests/test_integration_{steady,transient}.cpp`.
All mphys tests carry the ctest label `mphys`.

### Coverage

LLVM source-based coverage of `mphys_lib` (only `src/` and `include/mphys/` are
reported; external deps and test sources are excluded).

```bash
./scripts/coverage.sh            # build, run tests, print summary, open HTML report
./scripts/coverage.sh --no-open  # same, without opening the browser
```

Uses the `mac-coverage` CMake preset (`-fprofile-instr-generate
-fcoverage-mapping`). The build force-loads `mphys_lib` into the test binary
(`MPHYS_COVERAGE=ON`) so untested library code still counts against the total
instead of being dead-stripped by the linker. HTML report lands at
`build/mac-coverage/coverage/html/index.html`.

## Run

```bash
./build/mac-debug/gui/mphys_gui                          # Desktop GUI
./build/mac-debug/examples/example_convection_diffusion  # CLI example
```

## Code style

Google C++ Style Guide. All SUNDIALS objects wrapped in RAII structs (see `include/mphys/sundials_types.hpp`); never call SUNDIALS free functions manually.

## Architecture

**mphys_lib** is a cell-centred FVM library for 1D PDEs, built on SUNDIALS. Users subclass `Model`, declare fields, and implement `Residual()`. The same `Residual()` implementation works for both transient (IDA) and steady-state (KINSOL) solves — the distinction is whether `ydot` arrays are populated.

### Data flow for a simulation

```
Mesh1D          — geometry (cell centres, face positions, dx, df)
StateVector     — maps named Field objects to a flat SUNDIALS N_Vector (layout below)
Model           — holds fields[], algebraics[], bcs[]; user overrides Residual()
TransientSolver / SteadySolver — drives SUNDIALS, calls Residual() each step
SimResult       — sequence of Snapshots (t, fields[], algebraics[])
```

### StateVector memory layout

```
[var0_cell0, var1_cell0, ..., varK_cell0,
 var0_cell1, var1_cell1, ..., varK_cell1, ...
 alg0, alg1, ...]
```

Half-bandwidth of the banded Jacobian equals `n_vars`. Algebraic variables (0D scalars) are appended after all spatial data; they are marked with `id = 0` in the IDA id vector.

### FVM operators (`mphys::fvm::` namespace)

All operators return a new `Field`; they do not mutate their arguments.

| Operator | Notes |
|----------|-------|
| `Ddt(ydot)` | Identity marker — just returns `ydot` for readability in `Residual()` |
| `Grad(phi, mesh, bcs)` | Face-centred gradient, length `n_cells+1`; ghost cells from BCs |
| `Div(flux, mesh)` | Conservative divergence of a face-centred flux |
| `Laplacian(phi, D, mesh, bcs)` | Uniform or spatially-varying diffusivity (overloaded) |
| `Convection(phi, u, mesh, bcs)` | First-order upwind; uniform or face-varying velocity (overloaded) |
| `InterpolateToFaces(phi, mesh, bcs)` | Arithmetic-mean face interpolation |

### Boundary conditions

```cpp
SetBcs(field_idx, {mphys::DirichletBc(value), mphys::NeumannBc(value)});
// left BC = first argument, right BC = second argument
```

Dirichlet BCs are enforced via ghost cells inside the FVM operators; Neumann BCs set the face gradient directly.

### Adding a new physics model

1. Subclass `Model` in your `.cpp` or example file.
2. In the constructor: call `AddField("name", init_val)` for each spatial field, `AddAlgebraic("name", init_val)` for 0D scalars, then `SetBcs(idx, {left, right})`.
3. Implement `Residual()`: form the residual `rr[i] = ...` using `fvm::` operators. For steady-state, `ydot` arrays are empty — check before using `fvm::Ddt`.
4. Construct `StateVector sv(mesh.n_cells)`, instantiate your model, create `SolverOptions`, then run with `TransientSolver` or `SteadySolver`.

### GUI (`mphys_gui`)

COMSOL-inspired desktop app built with Dear ImGui (docking branch) + ImPlot + GLFW/OpenGL3. Requires `BUILD_GUI=ON` and `glfw` installed via Homebrew. Physics models are defined inline in `gui/mphys_gui.cpp`; the app stores results in `AppState::result` (a `SimResult`) alongside mesh `cell_centres` for plotting. The `imgui` submodule must use the `docking` branch — the master branch lacks the docking API.

## External dependencies

| Dependency | Location | Purpose |
|------------|----------|---------|
| SUNDIALS (IDA, IDAS, KINSOL) | `external/sundials` | ODE/DAE and nonlinear solvers |
| Google Test | `external/googletest` | Unit tests |
| Dear ImGui (docking branch) | `external/imgui` | GUI rendering |
| ImPlot | `external/implot` | Scientific plots in GUI |
| GLFW | system (`brew install glfw`) | Window + OpenGL context |
