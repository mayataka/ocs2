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
scalar_t QuadraticStoCost::getValue(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime, const ModeSchedule& modeSchedule, 
                                    const PreComputation&) const {
  const vector_t tsDeviation = getSwitchingTimeDeviation(initTime, switchingTimes, finalTime, modeSchedule);
  return 0.5 * tsDeviation.dot(Q_ * tsDeviation);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation QuadraticStoCost::getQuadraticApproximation(scalar_t initTime, const vector_t& switchingTimes, 
                                                                                 scalar_t finalTime, const ModeSchedule& modeSchedule, 
                                                                                 const PreComputation&) const {
  const vector_t tsDeviation = getSwitchingTimeDeviation(initTime, switchingTimes, finalTime, modeSchedule);

  ScalarFunctionQuadraticApproximation Phi;
  Phi.dfdxx = Q_;
  Phi.dfdx.noalias() = Q_ * tsDeviation;
  Phi.f = 0.5 * tsDeviation.dot(Phi.dfdx);
  return Phi;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t QuadraticStoCost::getSwitchingTimeDeviation(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime,  
                                                     const ModeSchedule& modeSchedule) const {
  vector_t tsDeviation(switchingTimes.size());
  for (size_t i = 0; i < switchingTimes.size(); ++i) {
    tsDeviation.coeffRef(i) = switchingTimes.coeff(i) - modeSchedule.eventTimes[i];
  }
  return tsDeviation;
}

}  // namespace ocs2
