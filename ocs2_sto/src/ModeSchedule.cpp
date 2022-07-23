#include <ocs2_sto/ModeSchedule.h>

namespace ocs2 {

scalar_array_t extractValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule) {
  const auto validModeSchedule = extractValidModeSchedule(initTime, finalTime, modeSchedule);
  return validModeSchedule.eventTimes;
}

ModeSchedule extractValidModeSchedule(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule) {
  ModeSchedule validModeSchedule;
  size_t lastIndex = 0;
  for (size_t i = 0 ; i < modeSchedule.eventTimes.size(); ++i) {
    lastIndex = i;
    if (modeSchedule.eventTimes[i] > initTime && modeSchedule.eventTimes[i] < finalTime) {
      validModeSchedule.eventTimes.push_back(modeSchedule.eventTimes[i]);
      validModeSchedule.modeSequence.push_back(modeSchedule.modeSequence[i]);
    }
    else if (modeSchedule.eventTimes[i] >= finalTime) break;
  }
  validModeSchedule.modeSequence.push_back(modeSchedule.modeSequence[lastIndex]);
  return validModeSchedule;
}

}  // namespace ocs2
