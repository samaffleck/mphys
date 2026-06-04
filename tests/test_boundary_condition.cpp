#include <gtest/gtest.h>

#include "mphys/boundary_condition.hpp"

namespace {

using mphys::BcType;
using mphys::BoundaryCondition;
using mphys::DirichletBc;
using mphys::FieldBcs;
using mphys::NeumannBc;

TEST(BoundaryCondition, DefaultIsZeroNeumann) {
  BoundaryCondition bc;
  EXPECT_EQ(bc.type, BcType::kNeumann);
  EXPECT_DOUBLE_EQ(bc.value, 0.0);
}

TEST(BoundaryCondition, DirichletFactory) {
  BoundaryCondition bc = DirichletBc(2.5);
  EXPECT_EQ(bc.type, BcType::kDirichlet);
  EXPECT_DOUBLE_EQ(bc.value, 2.5);
}

TEST(BoundaryCondition, NeumannFactory) {
  BoundaryCondition bc = NeumannBc(-1.0);
  EXPECT_EQ(bc.type, BcType::kNeumann);
  EXPECT_DOUBLE_EQ(bc.value, -1.0);
}

TEST(BoundaryCondition, FieldBcsHoldsBothSides) {
  FieldBcs bcs{DirichletBc(1.0), NeumannBc(3.0)};
  EXPECT_EQ(bcs.left.type, BcType::kDirichlet);
  EXPECT_DOUBLE_EQ(bcs.left.value, 1.0);
  EXPECT_EQ(bcs.right.type, BcType::kNeumann);
  EXPECT_DOUBLE_EQ(bcs.right.value, 3.0);
}

}  // namespace
