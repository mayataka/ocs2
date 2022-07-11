#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/core/SlackDual.h>
#include <ocs2_ipm/core/SlackDualDirection.h>

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

void updateIpmVariables(IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection,
                        scalar_t primalStepSize, scalar_t dualStepSize);

}  // namespace ipm
}  // namespace ocs2
