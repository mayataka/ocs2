#include <ocs2_sto/cost/QuadraticStoCost.h>

namespace ocs2 {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
QuadraticStoCost::QuadraticStoCost(matrix_t Q) : Q_(std::move(Q)), StoCost() {}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
QuadraticStoCost* QuadraticStoCost::clone() const {
  return new QuadraticStoCost(*this);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t QuadraticStoCost::getValue(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule, 
                                    const ModeSchedule& referenceModeSchedule, const PreComputation& preComp) const {
  const vector_t tsDeviation = getSwitchingTimeDeviation(initTime, finalTime, stoModeSchedule, referenceModeSchedule);
  return 0.5 * tsDeviation.dot(Q_ * tsDeviation);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation QuadraticStoCost::getQuadraticApproximation(scalar_t initTime, scalar_t finalTime, 
                                                                                 const ModeSchedule& stoModeSchedule, 
                                                                                 const ModeSchedule& referenceModeSchedule, 
                                                                                 const PreComputation& preComp) const {
  const vector_t tsDeviation = getSwitchingTimeDeviation(initTime, finalTime, stoModeSchedule, referenceModeSchedule);
  ScalarFunctionQuadraticApproximation Phi;
  Phi.dfdxx = Q_;
  Phi.dfdx.noalias() = Q_ * tsDeviation;
  Phi.f = 0.5 * tsDeviation.dot(Phi.dfdx);
  return Phi;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t QuadraticStoCost::getSwitchingTimeDeviation(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule, 
                                                     const ModeSchedule& referenceModeSchedule) const {
  const auto validSwitchingTimesPair = extractValidSwitchingTimes(initTime, finalTime, stoModeSchedule, referenceModeSchedule);
  const auto& validStoSwitchingTimes = validSwitchingTimesPair.first;
  const auto& validReferenceSwitchingTimes = validSwitchingTimesPair.second;
  vector_t tsDeviation(validStoSwitchingTimes.size());
  for (size_t i = 0; i < validStoSwitchingTimes.size(); ++i) {
    tsDeviation.coeffRef(i) = validStoSwitchingTimes[i] - validReferenceSwitchingTimes[i];
  }
  return tsDeviation;
}

}  // namespace ocs2
