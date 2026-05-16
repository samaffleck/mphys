/*
 * MPhys
 * Copyright (c) 2026 Sam Affleck and contributors
 *
 * This file is part of MPhys, an open source C++ library for multiphysics simulations.
 *
 * License: For academic and non-commercial use only.
 *
 * See the project README for more details.
 */
#include <cassert>
#include "mphys/solver_options.hpp"

namespace mphys {

  void SolverOptions::Check() const {
    assert(initial_time_step > 0.0);
    assert(minimum_time_step > 0.0);
    assert(maximum_time_step > 0.0);
    assert(minimum_time_step < initial_time_step);
    assert(minimum_time_step > 0.0);
    assert(maximum_time_step >= initial_time_step);
    assert(time_step_increase_factor >= 1.0);
    assert(time_step_decrease_factor >= 1.0);
    assert(jac_update_freq > 0);
  }
} // end mphys namespace
