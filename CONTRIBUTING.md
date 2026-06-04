# Contributing to mphys

Thanks for your interest in contributing! mphys is a C++23 finite volume
library for 1D PDEs, and contributions of all kinds are welcome — bug reports,
new physics models, FVM operators, documentation, and build/portability fixes.

## Getting started

1. Fork and clone the repository **with submodules** (SUNDIALS, ImGui, ImPlot,
   and Google Test are vendored as git submodules):

   ```bash
   git clone --recurse-submodules https://github.com/<you>/mphys.git
   # already cloned? grab the submodules:
   git submodule update --init --recursive
   ```

2. Configure and build (see the [README](README.md) for platform details):

   ```bash
   cmake --preset mac-debug        # or linux-debug
   cmake --build --preset mac-debug
   ```

3. Run the tests and the analytical validation suite before opening a PR:

   ```bash
   cd build/mac-debug && ctest --output-on-failure
   ./build/mac-debug/examples/example_validation_1d_diffusion
   ```

## Coding style

- Follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
- All SUNDIALS objects must be wrapped in the RAII types in
  `include/mphys/sundials_types.hpp` — never call SUNDIALS `*Free` functions
  directly.
- FVM operators in `mphys::fvm::` return a new `Field`; they never mutate their
  arguments.
- Keep public headers self-contained and documented with a short comment on
  each declaration, matching the existing style.

## Submitting changes

1. Create a topic branch off `master`.
2. Keep each PR focused on a single change; add or update tests where it makes
   sense.
3. Make sure the build is warning-clean and all tests pass.
4. Open a pull request with a clear description of *what* changed and *why*.
   Link any related issue.

## Adding a new physics model

See the "Defining a physics model" section of the [README](README.md) and the
worked examples in `examples/`. The recommended path is:

1. Subclass `mphys::Model`.
2. Declare fields/algebraics and boundary conditions in the constructor.
3. Implement `Residual()` using `mphys::fvm::` operators — the same
   implementation serves both transient (IDA) and steady-state (KINSOL) solves.

## Reporting bugs

Open an issue with a minimal reproduction: the governing equations, mesh and
boundary conditions, the observed vs. expected behaviour, and your platform and
compiler version.

## License

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE) that covers the project.
