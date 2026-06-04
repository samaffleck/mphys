# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

CMake presets are defined in `CMakePresets.json`, organised per platform — `mac`,
`linux`, `windows`, `web` — each with a `-debug` and `-release` variant (plus `-gui`
build presets for GUI-only iteration). Presets are guarded by `hostSystemName`, so only
the ones relevant to the current OS are offered. macOS uses Homebrew LLVM clang; the
other desktop presets use the system compiler. All desktop presets enable `BUILD_GUI=ON`.

```bash
# Configure + build (macOS shown; swap mac→linux/windows as appropriate)
cmake --preset mac-debug
cmake --preset mac-release
cmake --build --preset mac-debug
cmake --build --preset mac-release

# Build GUI only (fastest iteration on GUI changes)
cmake --build --preset mac-debug-gui
cmake --build --preset mac-release-gui
```

### WebAssembly build

The web build runs through Emscripten's `emcmake` wrapper (it injects
`CMAKE_TOOLCHAIN_FILE`, so the preset stays toolchain-agnostic and works with both
`emsdk` and Homebrew emscripten):

```bash
emcmake cmake --preset web-release
cmake --build --preset web-release-gui
# → build/web-release/gui/mphys_gui.{html,js,wasm,data}; serve over HTTP to run.
```

Cross-platform notes:
- **GLFW** is vendored as the `external/glfw` submodule and built from source unless a
  system GLFW is found (`find_package(glfw3 QUIET)` in the top-level `CMakeLists.txt`).
- The **web GUI** uses the Emscripten `contrib.glfw3` port + WebGL2 (GLES3); the render
  loop is driven by `emscripten_set_main_loop_arg`, and `gui/shell.html` is the HTML
  shell. Native file dialogs (tinyfiledialogs) and the assets path are guarded by
  `#ifdef __EMSCRIPTEN__` (assets are preloaded into the virtual FS at `/assets`).

VS Code launch configs ("mphys GUI (Debug/Release)", etc.) and matching build tasks are in `.vscode/`.

## Tests

```bash
# Build and run all tests
cmake --build build/mac-debug --target mphys_tests -j8
cd build/mac-debug && ctest --output-on-failure

# Run the validation suite (analytical convergence checks — more informative than gtest)
./build/mac-debug/examples/example_validation_1d_diffusion
```

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

COMSOL-inspired app built with Dear ImGui (docking branch) + ImPlot + GLFW/OpenGL3,
running on desktop (macOS/Windows/Linux) and in the browser via WebAssembly. Requires
`BUILD_GUI=ON`. Physics models are defined inline in `gui/mphys_gui.cpp`; the app stores
results in `AppState::result` (a `SimResult`) alongside mesh `cell_centres` for plotting.
The `imgui` submodule must use the `docking` branch — the master branch lacks the docking
API. The render loop is factored into `MainLoopStep()` so the same body serves both the
desktop `while` loop and the Emscripten `emscripten_set_main_loop_arg` callback.

## External dependencies

| Dependency | Location | Purpose |
|------------|----------|---------|
| SUNDIALS (IDA, IDAS, KINSOL) | `external/sundials` | ODE/DAE and nonlinear solvers |
| Google Test | `external/googletest` | Unit tests |
| Dear ImGui (docking branch) | `external/imgui` | GUI rendering |
| ImPlot | `external/implot` | Scientific plots in GUI |
| GLFW | `external/glfw` submodule (or system) | Window + OpenGL context (desktop) |
| tinyfiledialogs | `external/tinyfiledialogs` | Native file dialogs (desktop only) |

On the web build, GLFW + WebGL2 come from the Emscripten `contrib.glfw3` port instead of
the submodule; tinyfiledialogs is not compiled.
