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
scalar_t StoCostCollection::getValue(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime, const ModeSchedule& modeSchedule,
                                     const PreComputation& preComp) const {
  scalar_t cost = 0.0;

  // accumulate cost terms
  for (const auto& costTerm : this->terms_) {
    if (costTerm->isActive(initTime, switchingTimes, finalTime, modeSchedule)) {
      cost += costTerm->getValue(initTime, switchingTimes, finalTime, modeSchedule, preComp);
    }
  }

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation StoCostCollection::getQuadraticApproximation(scalar_t initTime, const vector_t& switchingTimes, 
                                                                                  scalar_t finalTime, const ModeSchedule& modeSchedule, 
                                                                                  const PreComputation& preComp) const {
  const auto firstActive =
      std::find_if(terms_.begin(), terms_.end(), [initTime, switchingTimes, finalTime, modeSchedule](const std::unique_ptr<StoCost>& costTerm) { 
          return costTerm->isActive(initTime, switchingTimes, finalTime, modeSchedule); 
      });

  // No active terms (or terms is empty).
  if (firstActive == terms_.end()) {
    return ScalarFunctionQuadraticApproximation::Zero(switchingTimes.size(), 0);
  }

  // Initialize with first active term, accumulate potentially other active terms.
  auto cost = (*firstActive)->getQuadraticApproximation(initTime, switchingTimes, finalTime, modeSchedule, preComp);
  std::for_each(std::next(firstActive), terms_.end(), [&](const std::unique_ptr<StoCost>& costTerm) {
    if (costTerm->isActive(initTime, switchingTimes, finalTime, modeSchedule)) {
      const auto costTermApproximation = costTerm->getQuadraticApproximation(initTime, switchingTimes, finalTime, modeSchedule, preComp);
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
