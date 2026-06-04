#include <gtest/gtest.h>

#include "mphys/solver_options.hpp"

namespace {

using mphys::SolverOptions;
using mphys::Tolerance;

TEST(SolverOptions, DefaultValues) {
  SolverOptions opts;
  EXPECT_DOUBLE_EQ(opts.initial_time_step, 1.0);
  EXPECT_DOUBLE_EQ(opts.minimum_time_step, 1e-10);
  EXPECT_DOUBLE_EQ(opts.maximum_time_step, 2.0);
  EXPECT_DOUBLE_EQ(opts.time_step_increase_factor, 2.0);
  EXPECT_DOUBLE_EQ(opts.time_step_decrease_factor, 2.0);
  EXPECT_DOUBLE_EQ(opts.relaxation_factor, 0.95);
  EXPECT_EQ(opts.jac_update_freq, 10);
  EXPECT_EQ(opts.max_Newton_iter, 200);
}

TEST(SolverOptions, DefaultToleranceValues) {
  Tolerance tol;
  EXPECT_DOUBLE_EQ(tol.absolute, 1e-12);
  EXPECT_DOUBLE_EQ(tol.relative, 1e-12);
}

TEST(SolverOptions, CheckAcceptsDefaultConfiguration) {
  SolverOptions opts;
  // Check() asserts on invalid configurations; the defaults must be valid.
  EXPECT_NO_FATAL_FAILURE(opts.Check());
}

TEST(SolverOptions, CheckAcceptsCustomValidConfiguration) {
  SolverOptions opts;
  opts.initial_time_step = 0.01;
  opts.minimum_time_step = 1e-8;
  opts.maximum_time_step = 0.5;
  opts.time_step_increase_factor = 1.5;
  opts.time_step_decrease_factor = 3.0;
  opts.jac_update_freq = 5;
  EXPECT_NO_FATAL_FAILURE(opts.Check());
}

}  // namespace
