/*
 * MPhys
 * Copyright (c) 2026 Sam Affleck and contributors
 *
 * This file is part of MPhys, an open source C++ library for multiphysics simulations.
 *
 * License: MIT. See the LICENSE file at the project root for details.
 */
#pragma once


namespace mphys {

  class Tolerance {
  public:
    double absolute = 1e-12;
    double relative = 1e-12;
  };

  class SolverOptions {
  public:
    void Check() const;

    Tolerance tolerance;
    double initial_time_step = 1.0;
    double minimum_time_step = 1e-10;
    double maximum_time_step = 2.0;
    double time_step_increase_factor = 2.0;
    double time_step_decrease_factor = 2.0;
    double relaxation_factor = 0.95;
    int jac_update_freq = 10;
    int max_Newton_iter = 200;
  };

} // end mphys namespace
