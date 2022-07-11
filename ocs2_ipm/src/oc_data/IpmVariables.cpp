#include <ocs2_ipm/oc_data/IpmVariables.h>

#include <cassert>

namespace ocs2 {
namespace ipm {

void updateIpmVariables(IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection,
                        scalar_t primalStepSize, scalar_t dualStepSize) {
  assert(primalStepSize >= 0.0);
  assert(primalStepSize <= 1.0);
  assert(dualStepSize >= 0.0);
  assert(dualStepSize <= 1.0);
  if (ipmVariables.slackDualStateIneqConstraint.slack.size() > 0) {
    assert(ipmVariables.slackDualStateIneqConstraint.slack.size() == ipmVariables.slackDualStateIneqConstraint.dual.size());
    ipmVariables.slackDualStateIneqConstraint.slack.noalias()
        += primalStepSize * ipmVariablesDirection.slackDualDirectionStateIneqConstraint.slackDirection; 
    ipmVariables.slackDualStateIneqConstraint.dual.noalias()
        += dualStepSize * ipmVariablesDirection.slackDualDirectionStateIneqConstraint.dualDirection; 
  }
  if (ipmVariables.slackDualStateInputIneqConstraint.slack.size() > 0) {
    assert(ipmVariables.slackDualStateInputIneqConstraint.slack.size() == ipmVariables.slackDualStateInputIneqConstraint.dual.size());
    ipmVariables.slackDualStateInputIneqConstraint.slack.noalias()
        += primalStepSize * ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint.slackDirection; 
    ipmVariables.slackDualStateInputIneqConstraint.dual.noalias()
        += dualStepSize * ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint.dualDirection; 
  }
}

}  // namespace ipm
}  // namespace ocs2
