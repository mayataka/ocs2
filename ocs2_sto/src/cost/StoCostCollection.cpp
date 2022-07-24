#include <ocs2_sto/cost/StoCostCollection.h>

namespace ocs2 {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
StoCostCollection::StoCostCollection(const StoCostCollection& other) : Collection<StoCost>(other) {}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
StoCostCollection* StoCostCollection::clone() const {
  return new StoCostCollection(*this);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t StoCostCollection::getValue(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule, 
                                     const ModeSchedule& referenceModeSchedule, const PreComputation& preComp) const {
  scalar_t cost = 0.0;

  // accumulate cost terms
  for (const auto& costTerm : this->terms_) {
    if (costTerm->isActive(initTime, finalTime, stoModeSchedule, referenceModeSchedule)) {
      cost += costTerm->getValue(initTime, finalTime, stoModeSchedule, referenceModeSchedule, preComp);
    }
  }

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation StoCostCollection::getQuadraticApproximation(scalar_t initTime, scalar_t finalTime, 
                                                                                  const ModeSchedule& stoModeSchedule, 
                                                                                  const ModeSchedule& referenceModeSchedule, 
                                                                                  const PreComputation& preComp) const {
  const auto firstActive =
      std::find_if(terms_.begin(), terms_.end(), [initTime, finalTime, stoModeSchedule, referenceModeSchedule](const std::unique_ptr<StoCost>& costTerm) { 
          return costTerm->isActive(initTime, finalTime, stoModeSchedule, referenceModeSchedule); 
      });

  // No active terms (or terms is empty).
  if (firstActive == terms_.end()) {
    return ScalarFunctionQuadraticApproximation::Zero(extractValidSwitchingTimes(initTime, finalTime, referenceModeSchedule).size(), 0);
  }

  // Initialize with first active term, accumulate potentially other active terms.
  auto cost = (*firstActive)->getQuadraticApproximation(initTime, finalTime, stoModeSchedule, referenceModeSchedule, preComp);
  std::for_each(std::next(firstActive), terms_.end(), [&](const std::unique_ptr<StoCost>& costTerm) {
    if (costTerm->isActive(initTime, finalTime, stoModeSchedule, referenceModeSchedule)) {
      const auto costTermApproximation = costTerm->getQuadraticApproximation(initTime, finalTime, stoModeSchedule, referenceModeSchedule, preComp);
      cost.f += costTermApproximation.f;
      cost.dfdx += costTermApproximation.dfdx;
      cost.dfdxx += costTermApproximation.dfdxx;
    }
  });

  // Make sure that input derivatives have zero size
  cost.dfdu.resize(0);
  cost.dfduu.resize(0, 0);
  cost.dfdux.resize(0, 0);

  return cost;
}

}  // namespace ocs2
