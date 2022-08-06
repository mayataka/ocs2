#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>
#include <ocs2_sto/model_data/StoModelData.h>
#include <ocs2_ipm/oc_problem/OptimalControlProblem.h>

namespace ocs2 {

/**
 * Calculates a quadratic approximate of the switching time optimization (STO) problem.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] initTime: The initial time.
 * @param [in] finalTime: The final time.
 * @param [in] referenceModeSchedule : The reference mode schedule.
 * @param [in] stoModeSchedule : The mode schedule for STO.
 * @param [out] stoModelData: The output data model.
 */
void approximateStoProblem(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                           const ModeSchedule& referenceModeSchedule, const ModeSchedule& stoModeSchedule, StoModelData& stoModelData);

/**
 * Calculates a quadratic approximate of the switching time optimization (STO) problem.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] initTime: The initial time.
 * @param [in] finalTime: The final time.
 * @param [in] referenceModeSchedule : The reference mode schedule.
 * @param [in] stoModeSchedule : The mode schedule for STO.
 * @return The output data model.
 */
StoModelData approximateStoProblem(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                                   const ModeSchedule& referenceModeSchedule, const ModeSchedule& stoModeSchedule);

/**
 * Compute the cost associated with the switching time optimization (STO). It is assumed that the precomputation request is already made.
 */
scalar_t computeCost(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                     const ModeSchedule& referenceModeSchedule, const ModeSchedule& stoModeSchedule);

/**
 * Compute the quadratic approximation of the cost associated with the switching time optimization (STO). It is assumed that the 
 * precomputation request is already made.
 */
ScalarFunctionQuadraticApproximation approximateCost(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                                                     const ModeSchedule& referenceModeSchedule, const ModeSchedule& stoModeSchedule);

}  // namespace ocs2
