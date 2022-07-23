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
size_t MinimumDwellTimeConstraint::getNumConstraints(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime, 
                                                     const ModeSchedule& /* modeSchedule */) const {
  const int firstActiveIndex = findFirstActiveIndex(initTime, switchingTimes, finalTime);
  const int lastActiveIndex = findLastActiveIndex(initTime, switchingTimes, finalTime);
  if (firstActiveIndex < 0 || lastActiveIndex < 0) {
    return 0;
  } else {
    assert(lastActiveIndex >= firstActiveIndex);
    return lastActiveIndex - firstActiveIndex + 2;
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t MinimumDwellTimeConstraint::getValue(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime, 
                                              const ModeSchedule& modeSchedule, const PreComputation& preComp) const {
  vector_t constraintValue(getNumConstraints(initTime, switchingTimes, finalTime, modeSchedule));
  const int firstActiveIndex = findFirstActiveIndex(initTime, switchingTimes, finalTime);
  const int lastActiveIndex = findLastActiveIndex(initTime, switchingTimes, finalTime);
  if (firstActiveIndex < 0 || lastActiveIndex < 0) {
    return constraintValue;
  } else {
    constraintValue.coeffRef(0) = switchingTimes.coeff(0) - 
    for (int i=firstActiveIndex; i <lastActiveIndex; ++i) {
      constraintValue.coeffRef(i)
    }
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximation MinimumDwellTimeConstraint::getLinearApproximation(scalar_t initTime, const vector_t& switchingTimes, 
                                                                                     scalar_t finalTime, const ModeSchedule& modeSchedule,
                                                                                     const PreComputation& preComp) const {
  const int firstActiveIndex = findFirstActiveIndex(initTime, switchingTimes, finalTime);
  const int lastActiveIndex = findLastActiveIndex(initTime, switchingTimes, finalTime);
  return 
}

} // namespace ocs2