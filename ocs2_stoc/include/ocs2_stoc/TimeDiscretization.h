#pragma once

#include <ocs2_core/NumericTraits.h>
#include <ocs2_core/Types.h>
#include <ocs2_sqp/TimeDiscretization.h>

namespace ocs2 {
namespace stoc {

/**
 * Decides on multi-phase time discretization along the horizon. Tries to makes step of dt, but will also ensure that eventtimes are part of the
 * discretization.
 *
 * @param initTime : start time.
 * @param finalTime : final time.
 * @param dt : desired discretization step.
 * @param eventTimes : Event times where a time discretization must be made.
 * @param dt_min : minimum discretization step. Smaller intervals will be merged. Needs to be bigger than limitEpsilon to avoid
 * interpolation problems
 * @return vector of discrete time points
 */
std::vector<AnnotatedTime> multiPhaseTimeDiscretization(scalar_t initTime, scalar_t finalTime, scalar_t dt,
                                                        const scalar_array_t& eventTimes,
                                                        scalar_t dt_min = 10.0 * numeric_traits::limitEpsilon<scalar_t>());

}  // namespace stoc 
}  // namespace ocs2
