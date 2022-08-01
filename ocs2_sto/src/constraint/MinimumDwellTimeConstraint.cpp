#include <ocs2_sto/constraint/MinimumDwellTimeConstraint.h>

#include <cassert>

namespace ocs2 {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
MinimumDwellTimeConstraint::MinimumDwellTimeConstraint(const std::unordered_map<size_t, scalar_t>& minimumDwellTimesMap,
                                                       const scalar_t minimumDwellTime) 
  : StoConstraint(), minimumDwellTimesMap_(minimumDwellTimesMap), minimumDwellTime_(minimumDwellTime) {}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
MinimumDwellTimeConstraint* MinimumDwellTimeConstraint::clone() const {
  return new MinimumDwellTimeConstraint(minimumDwellTimesMap_, minimumDwellTime_);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
size_t MinimumDwellTimeConstraint::getNumConstraints(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule,  
                                                     const ModeSchedule& /*stoModeSchedule*/) const { 
  const size_t numValidSwitchingTimes = getNumValidSwitchingTimes(initTime, finalTime, referenceModeSchedule);
  if (numValidSwitchingTimes == 0) {
    return 0;
  } else {
    return numValidSwitchingTimes + 1;
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t MinimumDwellTimeConstraint::getValue(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, 
                                              const ModeSchedule& stoModeSchedule, const PreComputation& preComp) const {
  const auto validModeSchedulePair = extractValidModeSchedulePair(initTime, finalTime, stoModeSchedule, referenceModeSchedule);
  const auto& validStoModeSchedule = validModeSchedulePair.first;
  if (validStoModeSchedule.eventTimes.empty()) {
    return vector_t();
  }
  const auto numConstraints = validStoModeSchedule.eventTimes.size() + 1;
  vector_t minimumDwellTimes(numConstraints);
  for (size_t i = 0; i < numConstraints; ++i) {
    const auto mode = validStoModeSchedule.modeSequence[i];
    const auto minDt = minimumDwellTimesMap_.find(mode);
    if (minDt != minimumDwellTimesMap_.end()) {
      minimumDwellTimes.coeffRef(i) = minDt->second;
    } else {
      minimumDwellTimes.coeffRef(i) = minimumDwellTime_;
    }
  }
  vector_t dwellTimes(numConstraints);
  dwellTimes.coeffRef(0) = validStoModeSchedule.eventTimes[0] - initTime;
  for (size_t i = 0; i < numConstraints - 2; ++i) {
    dwellTimes.coeffRef(i+1) = validStoModeSchedule.eventTimes[i+1] - validStoModeSchedule.eventTimes[i];
  }
  dwellTimes.coeffRef(numConstraints-1) = finalTime - validStoModeSchedule.eventTimes[numConstraints - 2];
  const vector_t constraintValue = minimumDwellTimes - dwellTimes;
  return constraintValue;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximation MinimumDwellTimeConstraint::getLinearApproximation(scalar_t initTime, scalar_t finalTime, 
                                                                                     const ModeSchedule& referenceModeSchedule, 
                                                                                     const ModeSchedule& stoModeSchedule, 
                                                                                     const PreComputation& preComp) const {
  VectorFunctionLinearApproximation linearApproximation(getNumConstraints(initTime, finalTime, referenceModeSchedule, stoModeSchedule), 
                                                        getNumValidSwitchingTimes(initTime, finalTime, referenceModeSchedule), 0);
  linearApproximation.f = getValue(initTime, finalTime, referenceModeSchedule, stoModeSchedule, preComp);
  const auto numConstraints = linearApproximation.f.size();
  if (numConstraints > 0) {
    linearApproximation.dfdx.setZero(numConstraints, numConstraints-1);
    linearApproximation.dfdx.coeffRef(0, 0) = -1.0;
    for (size_t i = 0; i < numConstraints - 2; ++i) {
      linearApproximation.dfdx.coeffRef(i+1, i)   =  1.0;
      linearApproximation.dfdx.coeffRef(i+1, i+1) = -1.0;
    }
    linearApproximation.dfdx.coeffRef(numConstraints-1, numConstraints-2) = 1.0;
  }
  return linearApproximation;
}

} // namespace ocs2