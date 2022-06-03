#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <ocs2_core/Types.h>
#include <ocs2_ipm/ipm/SlackDual.h>

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

}  // namespace ipm
}  // namespace ocs2
