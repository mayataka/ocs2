#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/core/SlackDual.h>
#include <ocs2_ipm/core/SlackDualDirection.h>
#include <ocs2_ipm/model_data/ModelData.h>

namespace ocs2 {
namespace ipm {

/**
 * The variables related to the interior point method.
 */
struct IpmVariables {
  SlackDual slackDualStateIneqConstraint; 
  SlackDual slackDualStateInputIneqConstraint; 
};

/**
 * The direction of the variables related to the interior point method.
 */
struct IpmVariablesDirection {
  SlackDualDirection slackDualDirectionStateIneqConstraint; 
  SlackDualDirection slackDualDirectionStateInputIneqConstraint; 
};

/**
 * Initializes the IPM variables (slack and dual variables).
 *
 * @param [in] modelData: The model data.
 * @param [out] ipmVariables: The ipm variables.
 * @param [in] barrier: The barrier parameter of the IPM.
 */
void initIpmVariables(const ModelData& modelData, IpmVariables& ipmVariables, scalar_t barrier);

/**
 * Initializes the IPM variables (slack and dual variables).
 *
 * @param [in] modelData: The model data.
 * @param [in] barrier: The barrier parameter of the IPM.
 * @return The ipm variables.
 */
IpmVariables initIpmVariables(const ModelData& modelData, scalar_t barrier);

/**
 * Initializes the IPM variables direction (directions of the slack and dual variables).
 *
 * @param [in] modelData: The model data.
 * @param [out] ipmVariablesDirection: The ipm variables direction.
 */
void initIpmVariablesDirection(const ModelData& modelData, IpmVariablesDirection& ipmVariablesDirection);

/**
 * Initializes the IPM variables direction (directions of the slack and dual variables).
 *
 * @param [in] ipmVariables: The ipm variables.
 * @param [out] ipmVariablesDirection: The ipm variables direction.
 */
void initIpmVariablesDirection(const IpmVariables& ipmVariables, IpmVariablesDirection& ipmVariablesDirection);

/**
 * Initializes the IPM variables direction (directions of the slack and dual variables).
 *
 * @param [in] modelData: The model data.
 * @return The ipm variables direction.
 */
IpmVariablesDirection initIpmVariablesDirection(const ModelData& modelData);

/**
 * Initializes the IPM variables direction (directions of the slack and dual variables).
 *
 * @param [in] ipmVariables: The ipm variables.
 * @return The ipm variables direction.
 */
IpmVariablesDirection initIpmVariablesDirection(const IpmVariables& ipmVariables);

/**
 * Updates the IPM variables iterate.
 *
 * @param [in, out] ipmVariables: The ipm variables.
 * @param [in] ipmVariablesDirection: The ipm variables direction.
 * @param [in] primalStepSize: The step size of the primal variables (slack variable).
 * @param [in] dualStepSize: The step size of the dual variables (slack variable).
 */
void updateIpmVariables(IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection,
                        scalar_t primalStepSize, scalar_t dualStepSize);

}  // namespace ipm
}  // namespace ocs2
