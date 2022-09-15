#include "ocs2_sto/constraint/MPC_MinimumDwellTimeConstraint.h"

#include <cassert>

namespace ocs2 {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
MPC_MinimumDwellTimeConstraint::MPC_MinimumDwellTimeConstraint(
    const scalar_t minimumDwellTimeAtInitialPhase, const scalar_t minimumDwellTimeAtFinalPhase,
    const std::unordered_map<size_t, scalar_t>& minimumDwellTimesMap, const scalar_t minimumDwellTime) 
  : MinimumDwellTimeConstraint(minimumDwellTimesMap, minimumDwellTime), 
    minimumDwellTimesMap_(minimumDwellTimesMap), 
    minimumDwellTime_(minimumDwellTime),
    minimumDwellTimeAtInitialPhase_(minimumDwellTimeAtInitialPhase),
    minimumDwellTimeAtFinalPhase_(minimumDwellTimeAtFinalPhase) {}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
MPC_MinimumDwellTimeConstraint* MPC_MinimumDwellTimeConstraint::clone() const {
  return new MPC_MinimumDwellTimeConstraint(minimumDwellTimeAtInitialPhase_, minimumDwellTimeAtFinalPhase_, minimumDwellTimesMap_, minimumDwellTime_);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t MPC_MinimumDwellTimeConstraint::getMinimumDwellTimes(const ModeSchedule& validModeSchedule) const {
  const auto numConstraints = validModeSchedule.eventTimes.size() + 1;
  vector_t minimumDwellTimes(numConstraints);
  for (size_t i = 0; i < numConstraints; ++i) {
    const auto mode = validModeSchedule.modeSequence[i];
    const auto minDt = minimumDwellTimesMap_.find(mode);
    if (minDt != minimumDwellTimesMap_.end()) {
      minimumDwellTimes.coeffRef(i) = minDt->second;
    } else {
      minimumDwellTimes.coeffRef(i) = minimumDwellTime_;
    }
  }
  minimumDwellTimes.coeffRef(0) = minimumDwellTimeAtInitialPhase_;
  minimumDwellTimes.coeffRef(minimumDwellTimes.size()-1) = minimumDwellTimeAtFinalPhase_;
  return minimumDwellTimes;
}

} // namespace ocs2