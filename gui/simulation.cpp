#include "simulation.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "mphys/field.hpp"
#include "mphys/fvm_mesh.hpp"
#include "mphys/mesh_model.hpp"
#include "mphys/mesh_steady_solver.hpp"
#include "mphys/mesh_transient_solver.hpp"
#include "mphys/models/model_registry.hpp"
#include "mphys/models/spm.hpp"
#include "mphys/models/spme.hpp"
#include "mphys/models/transport_mesh.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/steady_solver.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/transient_solver.hpp"

// ============================================================
// Geometry / parameter preparation
// ============================================================

// Returns cell-domain index for every cell in the composite 1D mesh.
static std::vector<int> CellDomainMap(const Geometry1D& geo) {
  std::vector<int> m;
  for (int d = 0; d < (int)geo.domains.size(); ++d)
    for (int i = 0; i < geo.domains[d].n_cells; ++i)
      m.push_back(d);
  return m;
}

const mphys::ModelInfo* EnsureModelConfig(AppState& s) {
  auto& reg = mphys::BuiltinModels();
  const mphys::ModelInfo* info = reg.Find(s.model_id);
  if (!info) {
    if (reg.All().empty()) return nullptr;
    info = &reg.All().front();
    s.model_id = info->id;
  }
  for (auto& dom : s.geo.domains)
    for (const auto& p : info->schema.params)
      if (p.scope == mphys::ParamScope::kPerDomain)
        dom.params.try_emplace(p.key, static_cast<float>(p.default_value));
  for (const auto& b : info->schema.boundaries)
    s.bcs.try_emplace(b.key, mphys::BcChoice{b.default_option, b.default_value});
  return info;
}

void ApplyLengthsToNodes(AppState& s) {
  s.geo.nodes.clear();
  s.geo.nodes.push_back({0.0f});
  float x = 0.0f;
  for (float L : s.geo_lengths) { x += L; s.geo.nodes.push_back({x}); }
  while ((int)s.geo.domains.size() < (int)s.geo_lengths.size())
    s.geo.domains.push_back({});
  s.geo.domains.resize(s.geo_lengths.size());
}

void ApplyNodesToLengths(AppState& s) {
  s.geo_lengths.clear();
  for (int d = 0; d < (int)s.geo.domains.size(); ++d)
    s.geo_lengths.push_back(s.geo.nodes[d + 1].x - s.geo.nodes[d].x);
}

// ============================================================
// Mesh builders
// ============================================================

static mphys::Mesh1D BuildCompositeMesh(const Geometry1D& geo, mphys::CoordSystem cs) {
  mphys::Mesh1D mesh;
  for (int d = 0; d < (int)geo.domains.size(); ++d) {
    double x0 = geo.nodes[d].x;
    double x1 = geo.nodes[d + 1].x;
    int    nc  = geo.domains[d].n_cells;
    auto   seg = mphys::MakeUniformMesh1D(x0, x1, nc, cs);
    for (int i = (d == 0 ? 0 : 1); i <= nc; ++i)
      mesh.face_positions.push_back(seg.face_positions[i]);
    for (auto v : seg.cell_centres) mesh.cell_centres.push_back(v);
    for (auto v : seg.dx)           mesh.dx.push_back(v);
  }
  mesh.n_cells = (int)mesh.cell_centres.size();
  mesh.df.resize(mesh.n_cells + 1);
  mesh.df[0] = 2.0 * (mesh.cell_centres[0] - mesh.face_positions[0]);
  for (int i = 1; i < mesh.n_cells; ++i)
    mesh.df[i] = mesh.cell_centres[i] - mesh.cell_centres[i - 1];
  mesh.df[mesh.n_cells] = 2.0 * (mesh.face_positions.back() - mesh.cell_centres.back());
  return mesh;
}

static mphys::Mesh BuildBoxMesh(const BoxGeometry& b) {
  const int nx = std::max(1, b.nx), ny = std::max(1, b.ny), nz = std::max(1, b.nz);
  if (b.dim == 3)
    return mphys::MakeStructuredMesh3D(b.x0, b.x1, nx, b.y0, b.y1, ny, b.z0, b.z1, nz);
  return mphys::MakeStructuredMesh2D(b.x0, b.x1, nx, b.y0, b.y1, ny);
}

// ============================================================
// Simulation runners
// ============================================================

// Single Particle Model: two normalised spherical particles + voltage output.
static void RunSpm(AppState& s) {
  try {
    auto p = ToSpmParameters(s.spm);
    int nc = std::max(2, s.spm.n_cells);

    auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, nc,
                                         mphys::CoordSystem::kSpherical);
    s.cell_centres = mesh.cell_centres;   // normalised radius r/R ∈ [0,1]

    mphys::StateVector sv(mesh.n_cells);
    mphys::models::SpmModel model(mesh, sv, p);

    mphys::SolverOptions opts;
    opts.tolerance.relative = static_cast<double>(s.rel_tol);
    opts.tolerance.absolute = static_cast<double>(s.abs_tol);
    opts.initial_time_step  = static_cast<double>(s.dt_initial);
    opts.maximum_time_step  = static_cast<double>(s.dt_max);

    double dt_snap = std::max(static_cast<double>(s.dt_snapshot),
                              static_cast<double>(s.dt_max));

    mphys::SunContext sunctx;
    mphys::TransientSolver solver(model, opts, sunctx);
    double next_snap = 0.0;
    std::string warn = solver.Solve(0.0, static_cast<double>(s.t_end),
        [&](double t, const std::vector<mphys::Field>& f,
            const std::vector<double>& a) {
          if (t >= next_snap - 1e-12) {
            s.result.Record(t, f, a);
            s.spm_time.push_back(t);
            s.spm_voltage.push_back(a[model.voltage_index()]);
            next_snap += dt_snap;
          }
        });

    if (!s.result.snapshots.empty()) {
      s.plot_time   = static_cast<float>(s.result.snapshots.back().t);
      s.has_results = true;
      s.nav         = NavNode::Results;
    }
    s.status_msg = warn.empty()
        ? "Done — " + std::to_string(s.result.snapshots.size()) + " snapshots"
        : "Warning: " + warn;
  } catch (const std::exception& e) {
    s.status_msg = std::string("Error: ") + e.what();
  }
}

// SPMe: two spherical particles + a macroscopic electrolyte field, plus the
// terminal-voltage output.
static void RunSpme(AppState& s) {
  try {
    auto p = ToSpmeParameters(s.spme);
    int nn = std::max(2, s.spme.n_n);
    int ns = std::max(1, s.spme.n_s);
    int np = std::max(2, s.spme.n_p);

    auto em = mphys::models::MakeSpmeElectrolyteMesh(p, nn, ns, np);
    auto particle_mesh = mphys::MakeUniformMesh1D(
        0.0, 1.0, em.mesh.n_cells, mphys::CoordSystem::kSpherical);

    s.cell_centres = particle_mesh.cell_centres;  // normalised radius r/R
    s.spme_ce_x    = em.mesh.cell_centres;         // electrolyte position x [m]

    mphys::StateVector sv(em.mesh.n_cells);
    mphys::models::SpmeModel model(particle_mesh, em, sv, p);

    mphys::SolverOptions opts;
    opts.tolerance.relative = static_cast<double>(s.rel_tol);
    opts.tolerance.absolute = static_cast<double>(s.abs_tol);
    opts.initial_time_step  = static_cast<double>(s.dt_initial);
    opts.maximum_time_step  = static_cast<double>(s.dt_max);

    double dt_snap = std::max(static_cast<double>(s.dt_snapshot),
                              static_cast<double>(s.dt_max));

    mphys::SunContext sunctx;
    mphys::TransientSolver solver(model, opts, sunctx);
    double next_snap = 0.0;
    std::string warn = solver.Solve(0.0, static_cast<double>(s.t_end),
        [&](double t, const std::vector<mphys::Field>& f,
            const std::vector<double>& a) {
          if (t >= next_snap - 1e-12) {
            s.result.Record(t, f, a);
            s.spm_time.push_back(t);
            s.spm_voltage.push_back(a[model.voltage_index()]);
            next_snap += dt_snap;
          }
        });

    if (!s.result.snapshots.empty()) {
      s.plot_time   = static_cast<float>(s.result.snapshots.back().t);
      s.has_results = true;
      s.nav         = NavNode::Results;
    }
    s.status_msg = warn.empty()
        ? "Done — " + std::to_string(s.result.snapshots.size()) + " snapshots"
        : "Warning: " + warn;
  } catch (const std::exception& e) {
    s.status_msg = std::string("Error: ") + e.what();
  }
}

// 2D/3D single-species convection-diffusion-reaction on a structured box,
// solved on the face-based mesh with the matrix-free solvers. Results are
// stored as per-cell field snapshots alongside the mesh for visualisation.
static void RunMeshSimulation(AppState& s) {
  const BoxGeometry& b = s.box;
  if (b.x1 <= b.x0 || b.y1 <= b.y0 || (b.dim == 3 && b.z1 <= b.z0))
    throw std::runtime_error("Box extents must be strictly increasing");

  s.mesh_result = BuildBoxMesh(b);

  const int npatch = static_cast<int>(s.mesh_result.patches.size());
  std::vector<mphys::PatchBc> bcs;
  bcs.reserve(npatch);
  for (int p = 0; p < npatch; ++p) {
    const int i = std::min(p, 5);
    bcs.push_back(b.bc_type[i] == 0 ? mphys::DirichletBc(b.bc_value[i])
                                    : mphys::NeumannBc(b.bc_value[i]));
  }

  mphys::models::ConvDiffReactionMesh model(
      s.mesh_result,
      {static_cast<double>(b.vx), static_cast<double>(b.vy),
       static_cast<double>(b.vz)},
      static_cast<double>(b.D), static_cast<double>(b.k), std::move(bcs));

  mphys::SolverOptions opts;
  opts.tolerance.relative = static_cast<double>(s.rel_tol);
  opts.tolerance.absolute = static_cast<double>(s.abs_tol);
  opts.initial_time_step  = static_cast<double>(s.dt_initial);
  opts.maximum_time_step  = static_cast<double>(s.dt_max);

  mphys::SunContext sunctx;
  std::string warn;
  if (b.steady) {
    mphys::MeshSteadySolver solver(model, opts, sunctx);
    solver.Solve();
    s.mesh_snaps.push_back({0.0, model.fields()[model.c()]});
  } else {
    mphys::MeshTransientSolver solver(model, opts, sunctx);
    double dt_snap = std::max(static_cast<double>(s.dt_snapshot),
                              static_cast<double>(s.dt_max));
    double next_snap = 0.0;
    warn = solver.Solve(0.0, static_cast<double>(s.t_end),
        [&](double t, const std::vector<std::vector<double>>& y) {
          if (t >= next_snap - 1e-12) {
            s.mesh_snaps.push_back({t, y[model.c()]});
            next_snap += dt_snap;
          }
        });
    // Guarantee a snapshot at the final time (model.fields() holds it).
    const double t_end = static_cast<double>(s.t_end);
    if (s.mesh_snaps.empty() || s.mesh_snaps.back().t < t_end - 1e-9)
      s.mesh_snaps.push_back({t_end, model.fields()[model.c()]});
  }

  s.result_dim = b.dim;
  if (!s.mesh_snaps.empty()) {
    s.plot_time   = static_cast<float>(s.mesh_snaps.back().t);
    s.has_results = true;
    s.nav         = NavNode::Results;
  }
  s.status_msg = warn.empty()
      ? "Done — " + std::to_string(s.mesh_snaps.size()) + " snapshot(s), " +
            std::to_string(s.mesh_result.NCells()) + " cells"
      : "Warning: " + warn;
}

void RunSimulation(AppState& s) {
  s.has_results = false;
  s.result.snapshots.clear();
  s.cell_centres.clear();
  s.spme_ce_x.clear();
  s.spm_time.clear();
  s.spm_voltage.clear();
  s.mesh_snaps.clear();
  s.result_dim = 1;
  s.status_msg = "Running...";

  try {
    const mphys::ModelInfo* info = EnsureModelConfig(s);
    if (!info) throw std::runtime_error("No physics model is registered");

    // Models that build their own geometry/results run through a custom path.
    if (s.model_id == kSpmId)  { RunSpm(s);  return; }
    if (s.model_id == kSpmeId) { RunSpme(s); return; }

    // 2D/3D convection-diffusion-reaction runs on the face-based mesh path.
    if (s.model_id == kConvDiffId && s.dim >= 2) { RunMeshSimulation(s); return; }

    if (!s.geo.built || s.geo.domains.empty())
      throw std::runtime_error("Build the geometry first (Geometry -> Build)");
    for (int d = 0; d < (int)s.geo.domains.size(); ++d)
      if (s.geo.nodes[d + 1].x <= s.geo.nodes[d].x)
        throw std::runtime_error("Node positions must be strictly increasing");

    auto mesh  = BuildCompositeMesh(s.geo, s.coord_system);
    s.cell_centres = mesh.cell_centres;
    auto cell_domain = CellDomainMap(s.geo);

    mphys::SolverOptions opts;
    opts.tolerance.relative = static_cast<double>(s.rel_tol);
    opts.tolerance.absolute = static_cast<double>(s.abs_tol);
    opts.initial_time_step  = static_cast<double>(s.dt_initial);
    opts.maximum_time_step  = static_cast<double>(s.dt_max);

    mphys::SunContext  sunctx;
    mphys::StateVector sv(mesh.n_cells);

    // The factory pulls per-domain parameters and boundary choices through
    // these lookups, expanding them to the per-cell / per-face data it needs.
    auto param_lookup = [&s](int d, const std::string& key) -> double {
      const auto& m = s.geo.domains[d].params;
      auto it = m.find(key);
      return it != m.end() ? static_cast<double>(it->second) : 0.0;
    };
    auto bc_lookup = [&s](const std::string& slot) -> mphys::BcChoice {
      auto it = s.bcs.find(slot);
      return it != s.bcs.end() ? it->second : mphys::BcChoice{};
    };
    mphys::ModelBuildContext ctx(mesh, sv, std::move(cell_domain),
                                 param_lookup, bc_lookup);
    auto model = info->factory(ctx);

    double dt_snap = std::max(static_cast<double>(s.dt_snapshot),
                               static_cast<double>(s.dt_max));
    std::string solver_warning;
    if (info->solver == mphys::SolverKind::kTransient) {
      mphys::TransientSolver solver(*model, opts, sunctx);
      double next_snap = 0.0;
      solver_warning = solver.Solve(0.0, static_cast<double>(s.t_end),
          [&](double t, const std::vector<mphys::Field>& f,
              const std::vector<double>& a) {
            if (t >= next_snap - 1e-12) {
              s.result.Record(t, f, a);
              next_snap += dt_snap;
            }
          });
    } else {
      mphys::SteadySolver solver(*model, opts, sunctx);
      solver.Solve();
      s.result.Record(0.0, model->fields(), model->algebraics());
    }

    if (!s.result.snapshots.empty()) {
      s.plot_time   = static_cast<float>(s.result.snapshots.back().t);
      s.has_results = true;
      s.nav         = NavNode::Results;
    }
    if (!solver_warning.empty())
      s.status_msg = "Warning: " + solver_warning +
                     " (" + std::to_string(s.result.snapshots.size()) + " snapshots)";
    else
      s.status_msg = "Done — " + std::to_string(s.result.snapshots.size()) + " snapshots";

  } catch (const std::exception& e) {
    s.status_msg = std::string("Error: ") + e.what();
  }
}
