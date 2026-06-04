#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include <nvector/nvector_serial.h>

#include "mphys/field.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/sun_context.hpp"

namespace {

using mphys::Field;
using mphys::StateVector;

TEST(StateVector, ReportsSizes) {
  StateVector sv(3);
  sv.AddField("a");
  sv.AddField("b");
  sv.AddAlgebraic("g");
  EXPECT_EQ(sv.NCells(), 3);
  EXPECT_EQ(sv.NVars(), 2);
  EXPECT_EQ(sv.NAlgebraics(), 1);
  EXPECT_EQ(sv.TotalSize(), 3 * 2 + 1);
}

TEST(StateVector, GatherUsesInterleavedLayout) {
  mphys::SunContext ctx;
  StateVector sv(3);
  sv.AddField("a");
  sv.AddField("b");
  sv.AddAlgebraic("g");

  mphys::SunVector nv(sv.TotalSize(), ctx);
  std::vector<Field> fields = {Field("a", 3), Field("b", 3)};
  fields[0].values = {10, 11, 12};
  fields[1].values = {20, 21, 22};
  std::vector<double> alg = {99};

  sv.Gather(fields, alg, nv);

  // Expected layout: [a0,b0, a1,b1, a2,b2, g]
  const double* d = N_VGetArrayPointer(nv.v);
  EXPECT_DOUBLE_EQ(d[0], 10);
  EXPECT_DOUBLE_EQ(d[1], 20);
  EXPECT_DOUBLE_EQ(d[2], 11);
  EXPECT_DOUBLE_EQ(d[3], 21);
  EXPECT_DOUBLE_EQ(d[4], 12);
  EXPECT_DOUBLE_EQ(d[5], 22);
  EXPECT_DOUBLE_EQ(d[6], 99);
}

TEST(StateVector, GatherScatterRoundTrip) {
  mphys::SunContext ctx;
  StateVector sv(4);
  sv.AddField("a");
  sv.AddField("b");
  sv.AddAlgebraic("g0");
  sv.AddAlgebraic("g1");

  std::vector<Field> in;
  std::vector<double> in_alg;
  sv.AllocateScratch(in, in_alg);
  for (int i = 0; i < 4; ++i) {
    in[0][i] = 1.0 * i;
    in[1][i] = 100.0 + i;
  }
  in_alg[0] = -7.0;
  in_alg[1] = 3.5;

  mphys::SunVector nv(sv.TotalSize(), ctx);
  sv.Gather(in, in_alg, nv);

  std::vector<Field> out;
  std::vector<double> out_alg;
  sv.AllocateScratch(out, out_alg);
  sv.Scatter(nv, out, out_alg);

  for (int k = 0; k < 2; ++k)
    for (int i = 0; i < 4; ++i) EXPECT_DOUBLE_EQ(out[k][i], in[k][i]);
  EXPECT_DOUBLE_EQ(out_alg[0], -7.0);
  EXPECT_DOUBLE_EQ(out_alg[1], 3.5);
}

TEST(StateVector, AllocateScratchSizesAndNamesMatch) {
  StateVector sv(5);
  sv.AddField("temperature");
  sv.AddField("pressure");
  sv.AddAlgebraic("flow");

  std::vector<Field> fields;
  std::vector<double> alg;
  sv.AllocateScratch(fields, alg);

  ASSERT_EQ(static_cast<int>(fields.size()), 2);
  EXPECT_EQ(fields[0].name, "temperature");
  EXPECT_EQ(fields[1].name, "pressure");
  EXPECT_EQ(fields[0].NCells(), 5);
  EXPECT_EQ(static_cast<int>(alg.size()), 1);
}

TEST(StateVector, FillIdVectorMarksDifferentialAndAlgebraic) {
  mphys::SunContext ctx;
  StateVector sv(2);
  sv.AddField("a");  // differential by default
  sv.AddField("b");
  sv.AddAlgebraic("g");
  sv.MarkFieldAlgebraic("b");  // make field b algebraic

  mphys::SunVector id(sv.TotalSize(), ctx);
  sv.FillIdVector(id);

  // Layout: [a0,b0, a1,b1, g] -> [1,0, 1,0, 0]
  const double* d = N_VGetArrayPointer(id.v);
  EXPECT_DOUBLE_EQ(d[0], 1.0);
  EXPECT_DOUBLE_EQ(d[1], 0.0);
  EXPECT_DOUBLE_EQ(d[2], 1.0);
  EXPECT_DOUBLE_EQ(d[3], 0.0);
  EXPECT_DOUBLE_EQ(d[4], 0.0);
}

TEST(StateVector, MarkFieldAlgebraicUnknownThrows) {
  StateVector sv(2);
  sv.AddField("a");
  EXPECT_THROW(sv.MarkFieldAlgebraic("does_not_exist"), std::invalid_argument);
}

}  // namespace
