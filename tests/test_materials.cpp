#include <algorithm>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "mphys/materials/database.hpp"
#include "mphys/materials/spm_materials.hpp"
#include "mphys/materials/spme_materials.hpp"
#include "mphys/mesh.hpp"
#include "mphys/models/spm.hpp"
#include "mphys/models/spme.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/transient_solver.hpp"

namespace {

using namespace mphys::materials;

template <typename T, typename U>
bool Contains(const std::vector<T>& v, const U& x) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

}  // namespace

// --- Registry / lookup ---

TEST(Materials, RegistryEnumeratesBuiltins) {
  EXPECT_GE(Database::ElectrodeIds().size(), 3u);
  EXPECT_GE(Database::ElectrolyteIds().size(), 1u);

  // Identifiers resolve to materials with the expected chemistry tags.
  EXPECT_EQ(Database::Electrode(ElectrodeId::kGraphiteChen2020).chemistry,
            "graphite");
  EXPECT_EQ(Database::Electrode(ElectrodeId::kNmc811Chen2020).chemistry,
            "NMC811");
  EXPECT_EQ(Database::Electrode(ElectrodeId::kLfpPrada2013).chemistry, "LFP");
}

TEST(Materials, DomainFilterSplitsElectrodes) {
  const auto neg = Database::ElectrodeIds(Domain::kNegativeElectrode);
  const auto pos = Database::ElectrodeIds(Domain::kPositiveElectrode);

  EXPECT_TRUE(Contains(neg, ElectrodeId::kGraphiteChen2020));
  EXPECT_FALSE(Contains(neg, ElectrodeId::kNmc811Chen2020));

  EXPECT_TRUE(Contains(pos, ElectrodeId::kNmc811Chen2020));
  EXPECT_TRUE(Contains(pos, ElectrodeId::kLfpPrada2013));
  EXPECT_FALSE(Contains(pos, ElectrodeId::kGraphiteChen2020));
}

TEST(Materials, DisplayNamesMatchIds) {
  EXPECT_EQ(Database::ElectrodeName(ElectrodeId::kNmc811Chen2020),
            "NMC811 (Chen 2020)");
  EXPECT_TRUE(Contains(Database::ElectrodeNames(), "Graphite (Chen 2020)"));
}

// --- Property sanity ---

TEST(Materials, ElectrodeOcpInPhysicalRange) {
  // Graphite sits low (< ~1.5 V), NMC/LFP high (> ~3 V) versus Li/Li+.
  const auto& graphite = Database::Electrode(ElectrodeId::kGraphiteChen2020);
  const auto& nmc = Database::Electrode(ElectrodeId::kNmc811Chen2020);
  const auto& lfp = Database::Electrode(ElectrodeId::kLfpPrada2013);

  EXPECT_GT(graphite.ocp(0.05), graphite.ocp(0.95));  // OCP falls on lithiation
  EXPECT_LT(graphite.ocp(0.5), 1.5);
  EXPECT_GT(nmc.ocp(0.5), 3.0);
  EXPECT_LT(nmc.ocp(0.5), 4.5);

  // LFP plateau: very flat between 20% and 80% lithiation.
  EXPECT_NEAR(lfp.ocp(0.2), lfp.ocp(0.8), 0.1);
}

TEST(Materials, ElectrolyteConductivityNearOneSiemensPerMetre) {
  const auto& e = Database::Electrolyte(ElectrolyteId::kLipf6EcEmcNyman2008);
  const double kappa = e.conductivity(1000.0, 298.15);   // 1 M
  EXPECT_GT(kappa, 0.5);
  EXPECT_LT(kappa, 1.5);
  const double De = e.diffusivity(1000.0, 298.15);
  EXPECT_GT(De, 1e-11);
  EXPECT_LT(De, 1e-9);
  EXPECT_NEAR(e.transference(1000.0, 298.15), 0.2594, 1e-9);
}

// --- Model integration: drive the SPM from a material selection ---

TEST(Materials, SpmRunsFromSelectedMaterials) {
  CellDesign cell;  // LG M50 defaults
  const mphys::models::SpmParameters p = MakeSpmParameters(
      ElectrodeId::kGraphiteChen2020, ElectrodeId::kNmc811Chen2020,
      ElectrolyteId::kLipf6EcEmcNyman2008, cell, /*current_A=*/5.0);

  // The adapter must map material data onto the parameter struct.
  EXPECT_DOUBLE_EQ(p.cn_max,
                   Database::Electrode(ElectrodeId::kGraphiteChen2020).c_max);
  EXPECT_DOUBLE_EQ(p.cp_max,
                   Database::Electrode(ElectrodeId::kNmc811Chen2020).c_max);
  EXPECT_DOUBLE_EQ(
      p.c_e, Database::Electrolyte(ElectrolyteId::kLipf6EcEmcNyman2008).c_typical);

  auto mesh = mphys::MakeUniformMesh1D(0.0, 1.0, 30,
                                       mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(mesh.n_cells);
  mphys::models::SpmModel model(mesh, sv, p);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-3;
  opts.maximum_time_step = 5.0;
  opts.tolerance.relative = 1e-8;
  opts.tolerance.absolute = 1e-8;

  mphys::SunContext sunctx;
  mphys::TransientSolver solver(model, opts, sunctx);

  std::vector<double> fa;
  solver.Solve(0.0, 50.0,
               [&](double, const std::vector<mphys::Field>&,
                   const std::vector<double>& a) { fa = a; });

  const double v = fa[model.voltage_index()];
  EXPECT_GT(v, 2.5);
  EXPECT_LT(v, 4.5);
}

TEST(Materials, SpmMaterialSlotsDescribeRequiredInputs) {
  const auto slots = SpmMaterialSlots();
  ASSERT_EQ(slots.size(), 3u);
  EXPECT_EQ(slots[0].category, Category::kElectrode);
  EXPECT_EQ(slots[0].domain, Domain::kNegativeElectrode);
  EXPECT_EQ(slots[2].category, Category::kElectrolyte);
}

// --- SPMe adapter ---

// The SPMe adapter must sample the electrolyte's concentration-dependent
// transport functions at the nominal salt concentration and copy them onto the
// constant-property SpmeParameters.
TEST(Materials, SpmeAdapterSamplesElectrolyteTransport) {
  const auto& elyte = Database::Electrolyte(ElectrolyteId::kLipf6EcEmcNyman2008);
  SpmeCellDesign cell;
  const auto p = MakeSpmeParameters(
      ElectrodeId::kGraphiteChen2020, ElectrodeId::kNmc811Chen2020,
      ElectrolyteId::kLipf6EcEmcNyman2008, cell, /*current_A=*/5.0);

  const double ce = elyte.c_typical;
  EXPECT_DOUBLE_EQ(p.ce0, ce);
  EXPECT_DOUBLE_EQ(p.D_e, elyte.diffusivity(ce, 298.15));
  EXPECT_DOUBLE_EQ(p.kappa_e, elyte.conductivity(ce, 298.15));
  EXPECT_DOUBLE_EQ(p.t_plus, elyte.transference(ce, 298.15));
  // Core particle data still flows through the SPM bridge.
  EXPECT_DOUBLE_EQ(p.core.cn_max,
                   Database::Electrode(ElectrodeId::kGraphiteChen2020).c_max);
}

TEST(Materials, SpmeRunsFromSelectedMaterials) {
  SpmeCellDesign cell;
  const auto p = MakeSpmeParameters(
      ElectrodeId::kGraphiteChen2020, ElectrodeId::kNmc811Chen2020,
      ElectrolyteId::kLipf6EcEmcNyman2008, cell, /*current_A=*/5.0);

  auto em = mphys::models::MakeSpmeElectrolyteMesh(p, 10, 5, 10);
  auto particle_mesh = mphys::MakeUniformMesh1D(
      0.0, 1.0, em.mesh.n_cells, mphys::CoordSystem::kSpherical);
  mphys::StateVector sv(em.mesh.n_cells);
  mphys::models::SpmeModel model(particle_mesh, em, sv, p);

  mphys::SolverOptions opts;
  opts.initial_time_step = 1e-3;
  opts.maximum_time_step = 5.0;
  opts.tolerance.relative = 1e-8;
  opts.tolerance.absolute = 1e-8;

  mphys::SunContext sunctx;
  mphys::TransientSolver solver(model, opts, sunctx);

  std::vector<double> fa;
  solver.Solve(0.0, 50.0,
               [&](double, const std::vector<mphys::Field>&,
                   const std::vector<double>& a) { fa = a; });

  const double v = fa[model.voltage_index()];
  EXPECT_GT(v, 2.5);
  EXPECT_LT(v, 4.5);
}
