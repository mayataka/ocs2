#include <ocs2_sto/ValidModeSchedule.h>

#include <cassert>

namespace ocs2 {

scalar_array_t extractValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule) {
  const auto validModeSchedule = extractValidModeSchedule(initTime, finalTime, modeSchedule);
  return std::move(validModeSchedule.eventTimes);
}

std::pair<scalar_array_t, scalar_array_t> extractValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, 
                                                                     const ModeSchedule& stoModeSchedule,
                                                                     const ModeSchedule& referenceModeSchedule) {
  ModeSchedule validStoModeSchedule, validReferenceModeSchedule;
  std::tie(validStoModeSchedule, validReferenceModeSchedule) = extractValidModeSchedule(initTime, finalTime, 
                                                                                        stoModeSchedule, referenceModeSchedule);
  return std::make_pair(std::move(validStoModeSchedule.eventTimes), std::move(validReferenceModeSchedule.eventTimes));
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

std::pair<ModeSchedule, ModeSchedule> extractValidModeSchedule(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule,
                                                               const ModeSchedule& referenceModeSchedule) {
  assert(stoModeSchedule.eventTimes.size() == referenceModeSchedule.eventTimes.size());
  assert(stoModeSchedule.modeSequence.size() == referenceModeSchedule.modeSequence.size());
  ModeSchedule validStoModeSchedule, validReferenceModeSchedule;
  size_t lastIndex = 0;
  for (size_t i = 0 ; i < referenceModeSchedule.eventTimes.size(); ++i) {
    lastIndex = i;
    if (referenceModeSchedule.eventTimes[i] > initTime && referenceModeSchedule.eventTimes[i] < finalTime) {
      validStoModeSchedule.eventTimes.push_back(stoModeSchedule.eventTimes[i]);
      validStoModeSchedule.modeSequence.push_back(stoModeSchedule.modeSequence[i]);

      validReferenceModeSchedule.eventTimes.push_back(referenceModeSchedule.eventTimes[i]);
      validReferenceModeSchedule.modeSequence.push_back(referenceModeSchedule.modeSequence[i]);
    }
    else if (referenceModeSchedule.eventTimes[i] >= finalTime) break;
  }
  validStoModeSchedule.modeSequence.push_back(stoModeSchedule.modeSequence[lastIndex]);
  validReferenceModeSchedule.modeSequence.push_back(referenceModeSchedule.modeSequence[lastIndex]);
  return std::make_pair(std::move(validStoModeSchedule), std::move(validReferenceModeSchedule));
}

}  // namespace ocs2
