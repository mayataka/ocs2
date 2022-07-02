#include <ocs2_ipm/oc_data/DualVariable.h>

namespace ocs2 {
namespace ipm {

void updateSlackDualIterate(DualVariable& dualVariable, const DualVariableDirection& dualDirection,
                            scalar_t primalStepSize, scalar_t dualStepSize) {
  if (dualVariable.slackDualStateIneqConstraint.slack.size() > 0) {
    dualVariable.slackDualStateIneqConstraint.slack.noalias()
        += primalStepSize * dualDirection.slackDualDirectionStateIneqConstraint.slackDirection; 
    dualVariable.slackDualStateIneqConstraint.dual.noalias()
        += dualStepSize * dualDirection.slackDualDirectionStateIneqConstraint.dualDirection; 
  }
  if (dualVariable.slackDualStateInputIneqConstraint.slack.size() > 0) {
    dualVariable.slackDualStateInputIneqConstraint.slack.noalias()
        += primalStepSize * dualDirection.slackDualDirectionStateInputIneqConstraint.slackDirection; 
    dualVariable.slackDualStateInputIneqConstraint.dual.noalias()
        += dualStepSize * dualDirection.slackDualDirectionStateInputIneqConstraint.dualDirection; 
  }
}

}  // namespace ipm
}  // namespace ocs2
