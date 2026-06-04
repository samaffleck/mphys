#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/state_vector.hpp"

namespace {

using mphys::BcType;
using mphys::DirichletBc;
using mphys::Field;
using mphys::FieldBcs;
using mphys::MakeUniformMesh1D;
using mphys::Mesh1D;
using mphys::Model;
using mphys::NeumannBc;
using mphys::StateVector;

// Minimal concrete model: Residual is required but irrelevant to these tests.
class TestModel : public Model {
 public:
  TestModel(const Mesh1D& mesh, StateVector& sv) : Model(mesh, sv) {}
  void Residual(double, const std::vector<Field>&, const std::vector<Field>&,
                const std::vector<double>&, const std::vector<double>&,
                std::vector<Field>& rr, std::vector<double>&) override {
    for (auto& r : rr)
      for (int i = 0; i < r.NCells(); ++i) r[i] = 0.0;
  }
};

class ModelTest : public ::testing::Test {
 protected:
  Mesh1D mesh_ = MakeUniformMesh1D(0.0, 1.0, 4);
  StateVector sv_{mesh_.n_cells};
  TestModel model_{mesh_, sv_};
};

TEST_F(ModelTest, AddFieldReturnsSequentialIndices) {
  EXPECT_EQ(model_.AddField("a"), 0);
  EXPECT_EQ(model_.AddField("b"), 1);
  EXPECT_EQ(model_.AddField("c"), 2);
  EXPECT_EQ(static_cast<int>(model_.fields().size()), 3);
}

TEST_F(ModelTest, AddFieldInitialisesValuesAndState) {
  model_.AddField("a", 5.0);
  EXPECT_EQ(model_.fields()[0].NCells(), mesh_.n_cells);
  for (int i = 0; i < mesh_.n_cells; ++i)
    EXPECT_DOUBLE_EQ(model_.fields()[0][i], 5.0);
  EXPECT_EQ(sv_.NVars(), 1);
}

TEST_F(ModelTest, NewFieldDefaultsToZeroNeumannBcs) {
  model_.AddField("a");
  const FieldBcs& bcs = model_.bcs()[0];
  EXPECT_EQ(bcs.left.type, BcType::kNeumann);
  EXPECT_EQ(bcs.right.type, BcType::kNeumann);
  EXPECT_DOUBLE_EQ(bcs.left.value, 0.0);
  EXPECT_DOUBLE_EQ(bcs.right.value, 0.0);
}

TEST_F(ModelTest, AddAlgebraicReturnsSequentialIndices) {
  EXPECT_EQ(model_.AddAlgebraic("p"), 0);
  EXPECT_EQ(model_.AddAlgebraic("q", 2.0), 1);
  EXPECT_EQ(static_cast<int>(model_.algebraics().size()), 2);
  EXPECT_DOUBLE_EQ(model_.algebraics()[1], 2.0);
  EXPECT_EQ(sv_.NAlgebraics(), 2);
}

TEST_F(ModelTest, SetBcsByIndex) {
  model_.AddField("a");
  model_.SetBcs(0, {DirichletBc(1.0), NeumannBc(2.0)});
  EXPECT_EQ(model_.bcs()[0].left.type, BcType::kDirichlet);
  EXPECT_DOUBLE_EQ(model_.bcs()[0].left.value, 1.0);
  EXPECT_DOUBLE_EQ(model_.bcs()[0].right.value, 2.0);
}

TEST_F(ModelTest, SetBcsByName) {
  model_.AddField("a");
  model_.AddField("b");
  model_.SetBcs("b", {DirichletBc(9.0), DirichletBc(8.0)});
  EXPECT_DOUBLE_EQ(model_.bcs()[1].left.value, 9.0);
  EXPECT_DOUBLE_EQ(model_.bcs()[1].right.value, 8.0);
  // Field "a" untouched.
  EXPECT_EQ(model_.bcs()[0].left.type, BcType::kNeumann);
}

TEST_F(ModelTest, SetBcsByIndexOutOfRangeThrows) {
  model_.AddField("a");
  EXPECT_THROW(model_.SetBcs(5, {DirichletBc(0.0), DirichletBc(0.0)}),
               std::out_of_range);
  EXPECT_THROW(model_.SetBcs(-1, {DirichletBc(0.0), DirichletBc(0.0)}),
               std::out_of_range);
}

TEST_F(ModelTest, SetBcsByUnknownNameThrows) {
  model_.AddField("a");
  EXPECT_THROW(model_.SetBcs("missing", {DirichletBc(0.0), DirichletBc(0.0)}),
               std::invalid_argument);
}

TEST_F(ModelTest, MeshAccessorReturnsConstructionMesh) {
  EXPECT_EQ(model_.mesh().n_cells, mesh_.n_cells);
}

}  // namespace
