#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/ipm/SlackDual.h>
#include <ocs2_ipm/ipm/SlackDualDirection.h>

namespace ocs2 {
namespace ipm {

/**
 * The dual variables for the interior point method.
 */
struct DualVariable {
  // Costate, i.e., the Lagrange multiplier w.r.t. the discretized state equation
  vector_t costate;

  // Interior point method related variables
  SlackDual slackDualStateIneqConstraint; 
  SlackDual slackDualStateInputIneqConstraint; 
};

/**
 * The direction of the dual variables for the interior point method.
 */
struct DualVariableDirection {
  // Interior point method related variables
  SlackDualDirection slackDualDirectionStateIneqConstraint; 
  SlackDualDirection slackDualDirectionStateInputIneqConstraint; 
};

}  // namespace ipm
}  // namespace ocs2
