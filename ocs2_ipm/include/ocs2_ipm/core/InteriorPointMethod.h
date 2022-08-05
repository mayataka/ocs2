#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/core/SlackDual.h>
#include <ocs2_ipm/core/SlackDualDirection.h>
#include <ocs2_ipm/core/InteriorPointMethodData.h>

namespace ocs2 {
namespace ipm {

/**
 * Initializes the slack and dual variables.
 *
 * @param [in] ineqConstraint: The inequality constraint.
 * @param [out] slackDual: Slack and dual variables of the inequality constraint.
 * @param [in] barrierParam: The barrier parameter of the IPM.
 */
void initSlackDual(const VectorFunctionLinearApproximation& ineqConstraint, SlackDual& slackDual, scalar_t barrierParam);

/**
 * Initializes the slack and dual variables.
 *
 * @param [in] ineqConstraint: The inequality constraint.
 * @param [in] barrierParam: The barrier parameter of the IPM.
 * @return Slack and dual variables of the inequality constraint.
 */
SlackDual initSlackDual(const VectorFunctionLinearApproximation& ineqConstraint, scalar_t barrierParam);

/**
 * Initializes the IPM data.
 *
 * @param [in] ineqConstraint: The inequality constraint.
 * @param [out] ipmData: IPM data of the inequality constraint.
 */
void initInteriorPointMethodData(const VectorFunctionLinearApproximation& ineqConstraint, InteriorPointMethodData& ipmData);

/**
 * Initializes the IPM data.
 *
 * @param [in] ineqConstraint: The inequality constraint.
 * @return IPM data of the inequality constraint.
 */
InteriorPointMethodData initInteriorPointMethodData(const VectorFunctionLinearApproximation& ineqConstraint);

/**
 * Evaluates the perturbed IPM residuals, i.e., the primal residual, perturbed complementary slackness, and contribution to the
 * derivatives of the Lagrangian with respect to the primal variables. 
 *
 * @param [in] ineqConstraint: The inequality constraint.
 * @param [in] slackDual: The slack and dual variables of the inequality constraint.
 * @param [out] ipmData: IPM data of the inequality constraint.
 * @param [in, out] cost: The quadratic model of the derivatives of the Lagrangian with respect to the primal variables.
 * @param [in] barrierParam: The barrier parameter of the IPM.
 * @param [in] stateConstraint: Flaf if the state-only inequality constraint is active.
 * @param [in] inputConstraint: Flaf if the state-input inequality constraint is active.
 */
void evalPerturbedResidual(const VectorFunctionLinearApproximation& ineqConstraint, const SlackDual& slackDual, 
                           InteriorPointMethodData& ipmData, ScalarFunctionQuadraticApproximation& cost, scalar_t barrierParam,
                           bool stateConstraint=true, bool inputConstraint=true);

/**
 * Condenses the inequality constraints and IPM-related variables from the quadratic model of the Lagrangian. It is assumed that the
 * perturbed IPM residuals are already evaluated.
 *
 * @param [in] ineqConstraint: The inequality constraint.
 * @param [in] slackDual: The slack and dual variables of the inequality constraint.
 * @param [out] ipmData: IPM data of the inequality constraint.
 * @param [in, out] cost: The quadratic model of the derivatives of the Lagrangian with respect to the primal variables.
 * @param [in] stateConstraint: Flaf if the state-only inequality constraint is active.
 * @param [in] inputConstraint: Flaf if the state-input inequality constraint is active.
 */
void condenseSlackDual(const VectorFunctionLinearApproximation& ineqConstraint, const SlackDual& slackDual, 
                       InteriorPointMethodData& ipmData, ScalarFunctionQuadraticApproximation& cost,
                       bool stateConstraint=true, bool inputConstraint=true);

/**
 * Expands the direction of the dual variable of the IPM method. It is assumed that the slack direction is already computed.
 *
 * @param [in] slackDual: The slack and dual variables of the inequality constraint.
 * @param [in] ipmData: IPM data of the inequality constraint.
 * @param [in, out] slackDualDirection: Directions of the slack and dual variables.
 */
void expandDual(const SlackDual& slackDual, const InteriorPointMethodData& ipmData, 
                SlackDualDirection& slackDualDirection);

/**
 * Expands the direction of the slack and dual variable of the IPM method from the directions of th primal variables.
 *
 * @param [in] ineqConstraint: The inequality constraint.
 * @param [in] slackDual: The slack and dual variables of the inequality constraint.
 * @param [in] ipmData: IPM data of the inequality constraint.
 * @param [in] dx: The state direction.
 * @param [in] du: The input direction.
 * @param [in, out] slackDualDirection: Directions of the slack and dual variables.
 */
void expandSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                     const SlackDual& slackDual, const InteriorPointMethodData& ipmData, 
                     const vector_t& dx, const vector_t& du, SlackDualDirection& slackDualDirection);

/**
 * Expands the direction of the slack and dual variable of the IPM method from the directions of th primal variables.
 *
 * @param [in] ineqConstraint: The inequality constraint.
 * @param [in] slackDual: The slack and dual variables of the inequality constraint.
 * @param [in] ipmData: IPM data of the inequality constraint.
 * @param [in] dx: The state direction.
 * @param [in, out] slackDualDirection: Directions of the slack and dual variables.
 */
void expandSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                     const SlackDual& slackDual, const InteriorPointMethodData& ipmData, 
                     const vector_t& dx, SlackDualDirection& slackDualDirection);

/**
 * Computes the step size via the fraction-to-boundary rule.
 *
 * @param [in] dim: Dimension of the vector.
 * @param [in] v: The vector.
 * @param [in] dv: The direction of the vector.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The step size.
 */
scalar_t fractionToBoundaryStepSize(size_t dim, const vector_t& v, const vector_t& dv,
                                    scalar_t marginRate=0.995);

/**
 * Computes the primal step size via the fraction-to-boundary rule.
 *
 * @param [in] slackDual: The slack and dual variables of the inequality constraint.
 * @param [in] slackDualDirection: Directions of the slack and dual variables.
 * @param [in] ipmData: IPM data of the inequality constraint.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The primal step size.
 */
scalar_t fractionToBoundaryPrimalStepSize(const SlackDual& slackDual, const SlackDualDirection& slackDualDirection,
                                          const InteriorPointMethodData& ipmData, scalar_t marginRate=0.995);

/**
 * Computes the dual step size via the fraction-to-boundary rule.
 *
 * @param [in] slackDual: The slack and dual variables of the inequality constraint.
 * @param [in] slackDualDirection: Directions of the slack and dual variables.
 * @param [in] ipmData: IPM data of the inequality constraint.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The dual step size.
 */
scalar_t fractionToBoundaryDualStepSize(const SlackDual& slackDual, const SlackDualDirection& slackDualDirection,
                                        const InteriorPointMethodData& ipmData, scalar_t marginRate=0.995);

}  // namespace ipm
}  // namespace ocs2
