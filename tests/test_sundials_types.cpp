// Unit tests for the RAII SUNDIALS wrappers. These guarantee single ownership
// and safe move semantics so that solver objects never double-free.
#include <gtest/gtest.h>

#include <utility>

#include "mphys/sundials_types.hpp"

namespace {

using mphys::IdaMem;
using mphys::KinMem;
using mphys::SunContext;
using mphys::SunLinearSolver;
using mphys::SunMatrix;
using mphys::SunVector;

TEST(SundialsTypes, ContextCreatesHandle) {
  SunContext ctx;
  EXPECT_NE(ctx.ctx, nullptr);
  EXPECT_NE(static_cast<SUNContext>(ctx), nullptr);
}

TEST(SundialsTypes, ContextMoveTransfersOwnership) {
  SunContext a;
  SUNContext raw = a.ctx;
  SunContext b(std::move(a));
  EXPECT_EQ(b.ctx, raw);
  EXPECT_EQ(a.ctx, nullptr);  // moved-from is null, safe to destroy

  SunContext c;
  c = std::move(b);
  EXPECT_EQ(c.ctx, raw);
  EXPECT_EQ(b.ctx, nullptr);
}

TEST(SundialsTypes, VectorAllocatesAndExposesData) {
  SunContext ctx;
  SunVector v(4, ctx);
  ASSERT_NE(v.v, nullptr);
  v.Data()[0] = 1.5;
  v.Data()[3] = -2.0;
  const SunVector& cv = v;
  EXPECT_DOUBLE_EQ(cv.Data()[0], 1.5);
  EXPECT_DOUBLE_EQ(cv.Data()[3], -2.0);
}

TEST(SundialsTypes, VectorMoveTransfersOwnership) {
  SunContext ctx;
  SunVector a(3, ctx);
  N_Vector raw = a.v;
  SunVector b(std::move(a));
  EXPECT_EQ(b.v, raw);
  EXPECT_EQ(a.v, nullptr);

  SunVector c;
  c = std::move(b);
  EXPECT_EQ(c.v, raw);
  EXPECT_EQ(b.v, nullptr);
}

TEST(SundialsTypes, MatrixMoveTransfersOwnership) {
  SunContext ctx;
  SunMatrix a(5, 2, 2, ctx);
  SUNMatrix raw = a.m;
  SunMatrix b(std::move(a));
  EXPECT_EQ(b.m, raw);
  EXPECT_EQ(a.m, nullptr);

  SunMatrix c;
  c = std::move(b);
  EXPECT_EQ(c.m, raw);
  EXPECT_EQ(b.m, nullptr);
}

TEST(SundialsTypes, LinearSolverMoveTransfersOwnership) {
  SunContext ctx;
  SunVector v(5, ctx);
  SunMatrix m(5, 2, 2, ctx);
  SunLinearSolver a(v, m, ctx);
  SUNLinearSolver raw = a.ls;
  SunLinearSolver b(std::move(a));
  EXPECT_EQ(b.ls, raw);
  EXPECT_EQ(a.ls, nullptr);

  SunLinearSolver c;
  c = std::move(b);
  EXPECT_EQ(c.ls, raw);
  EXPECT_EQ(b.ls, nullptr);
}

TEST(SundialsTypes, IdaMemMoveTransfersOwnership) {
  SunContext ctx;
  IdaMem a(ctx);
  void* raw = a.mem;
  IdaMem b(std::move(a));
  EXPECT_EQ(b.mem, raw);
  EXPECT_EQ(a.mem, nullptr);

  IdaMem c;
  c = std::move(b);
  EXPECT_EQ(c.mem, raw);
  EXPECT_EQ(b.mem, nullptr);
}

TEST(SundialsTypes, KinMemMoveTransfersOwnership) {
  SunContext ctx;
  KinMem a(ctx);
  void* raw = a.mem;
  KinMem b(std::move(a));
  EXPECT_EQ(b.mem, raw);
  EXPECT_EQ(a.mem, nullptr);

  KinMem c;
  c = std::move(b);
  EXPECT_EQ(c.mem, raw);
  EXPECT_EQ(b.mem, nullptr);
}

}  // namespace
