#pragma once

// Umbrella header: include this to pull in the entire public mphys API.
//
//   #include <mphys/mphys.hpp>
//
// Individual headers may also be included directly if you prefer finer-grained
// dependencies.

#include "mphys/boundary_condition.hpp"
#include "mphys/field.hpp"
#include "mphys/fvm_operators.hpp"
#include "mphys/mesh.hpp"
#include "mphys/model.hpp"
#include "mphys/sim_result.hpp"
#include "mphys/solver_options.hpp"
#include "mphys/state_vector.hpp"
#include "mphys/steady_solver.hpp"
#include "mphys/sun_context.hpp"
#include "mphys/sundials_types.hpp"
#include "mphys/transient_solver.hpp"

namespace mphys {

inline constexpr int version_major = 1;
inline constexpr int version_minor = 0;

}  // namespace mphys
