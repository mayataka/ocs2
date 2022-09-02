#include "ocs2_ipm_oc/oc_data/IpmVariables.h"

#include <ocs2_ipm/core/InteriorPointMethod.h>

#include <cassert>

namespace ocs2 {
namespace ipm {

void initIpmVariables(const ModelData& modelData, IpmVariables& ipmVariables, scalar_t barrier) {
  initSlackDual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, barrier);
  initSlackDual(modelData.stateInputIneqConstraint, ipmVariables.slackDualStateInputIneqConstraint, barrier);
}

IpmVariables initIpmVariables(const ModelData& modelData, scalar_t barrier) {
  IpmVariables ipmVariables;
  initIpmVariables(modelData, ipmVariables, barrier);
  return ipmVariables;
}

void initIpmVariablesDirection(const ModelData& modelData, IpmVariablesDirection& ipmVariablesDirection) {
  ipmVariablesDirection.slackDualDirectionStateIneqConstraint.resize(modelData.stateIneqConstraint.f.size());
  ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint.resize(modelData.stateInputIneqConstraint.f.size());
}

void initIpmVariablesDirection(const IpmVariables& ipmVariables, IpmVariablesDirection& ipmVariablesDirection) {
  ipmVariablesDirection.slackDualDirectionStateIneqConstraint.resize(ipmVariables.slackDualStateIneqConstraint.slack.size());
  ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint.resize(ipmVariables.slackDualStateInputIneqConstraint.slack.size());
}

IpmVariablesDirection initIpmVariablesDirection(const ModelData& modelData) {
  IpmVariablesDirection ipmVariablesDirection;
  initIpmVariablesDirection(modelData, ipmVariablesDirection);
  return ipmVariablesDirection;
}

IpmVariablesDirection initIpmVariablesDirection(const IpmVariables& ipmVariables) {
  IpmVariablesDirection ipmVariablesDirection;
  initIpmVariablesDirection(ipmVariables, ipmVariablesDirection);
  return ipmVariablesDirection;
}

void updateIpmVariables(IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection,
                        scalar_t primalStepSize, scalar_t dualStepSize) {
  assert(primalStepSize >= 0.0);
  assert(primalStepSize <= 1.0);
  assert(dualStepSize >= 0.0);
  assert(dualStepSize <= 1.0);
  if (ipmVariables.slackDualStateIneqConstraint.slack.size() > 0) {
    assert(ipmVariables.slackDualStateIneqConstraint.slack.size() == ipmVariables.slackDualStateIneqConstraint.dual.size());
    assert(ipmVariables.slackDualStateIneqConstraint.slack.size() 
            == ipmVariablesDirection.slackDualDirectionStateIneqConstraint.slackDirection.size());
    assert(ipmVariables.slackDualStateIneqConstraint.dual.size() 
            == ipmVariablesDirection.slackDualDirectionStateIneqConstraint.dualDirection.size());
    ipmVariables.slackDualStateIneqConstraint.slack.noalias()
        += primalStepSize * ipmVariablesDirection.slackDualDirectionStateIneqConstraint.slackDirection; 
    ipmVariables.slackDualStateIneqConstraint.dual.noalias()
        += dualStepSize * ipmVariablesDirection.slackDualDirectionStateIneqConstraint.dualDirection; 
  }
  if (ipmVariables.slackDualStateInputIneqConstraint.slack.size() > 0) {
    assert(ipmVariables.slackDualStateInputIneqConstraint.slack.size() == ipmVariables.slackDualStateInputIneqConstraint.dual.size());
    assert(ipmVariables.slackDualStateInputIneqConstraint.slack.size() 
            == ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint.slackDirection.size());
    assert(ipmVariables.slackDualStateInputIneqConstraint.dual.size() 
            == ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint.dualDirection.size());
    ipmVariables.slackDualStateInputIneqConstraint.slack.noalias()
        += primalStepSize * ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint.slackDirection; 
    ipmVariables.slackDualStateInputIneqConstraint.dual.noalias()
        += dualStepSize * ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint.dualDirection; 
  }
}

}  // namespace ipm
}  // namespace ocs2
