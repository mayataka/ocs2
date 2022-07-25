#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>
#include <ocs2_sto/model_data/StoModelData.h>
#include <ocs2_ipm/oc_problem/OptimalControlProblem.h>

namespace ocs2 {
namespace stoc {

void approximateStoProblem(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                           const ModeSchedule& stoModeSchedule, const ModeSchedule& referenceModeSchedule, StoModelData& stoModelData);

StoModelData approximateStoProblem(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                                   const ModeSchedule& stoModeSchedule, const ModeSchedule& referenceModeSchedule);

/**
 * Compute the total intermediate cost (i.e. cost + softConstraints). It is assumed that the precomputation request is already made.
 */
scalar_t computeCost(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                     const ModeSchedule& stoModeSchedule, const ModeSchedule& referenceModeSchedule);

/**
 * Compute the quadratic approximation of the total intermediate cost (i.e. cost + softConstraints). It is assumed that the precomputation
 * request is already made.
 */
ScalarFunctionQuadraticApproximation approximateCost(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                                                     const ModeSchedule& stoModeSchedule, const ModeSchedule& referenceModeSchedule);

}  // namespace stoc
}  // namespace ocs2
