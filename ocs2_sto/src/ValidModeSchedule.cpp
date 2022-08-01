#include <ocs2_sto/ValidModeSchedule.h>

#include <cassert>

namespace ocs2 {

size_t getNumValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule) {
  size_t num = 0;
  for (size_t i = 0 ; i < modeSchedule.eventTimes.size(); ++i) {
    num = i;
    if (modeSchedule.eventTimes[i] >= finalTime) break;
  }
  return num;
}

size_array_t extractValidSwitchingTimeIndices(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule) {
  size_array_t validSwitchingTimeIndices;
  size_t lastIndex = 0;
  for (size_t i = 0 ; i < modeSchedule.eventTimes.size(); ++i) {
    lastIndex = i;
    if (modeSchedule.eventTimes[i] > initTime && modeSchedule.eventTimes[i] < finalTime) {
      validSwitchingTimeIndices.push_back(i);
    }
    else if (modeSchedule.eventTimes[i] >= finalTime) break;
  }
  return validSwitchingTimeIndices;
}

scalar_array_t extractValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule) {
  const auto validModeSchedule = extractValidModeSchedule(initTime, finalTime, modeSchedule);
  return std::move(validModeSchedule.eventTimes);
}

std::pair<scalar_array_t, scalar_array_t> extractValidSwitchingTimesPair(scalar_t initTime, scalar_t finalTime, 
                                                                         const ModeSchedule& referenceModeSchedule,
                                                                         const ModeSchedule& stoModeSchedule) {
  ModeSchedule validReferenceModeSchedule, validStoModeSchedule;
  std::tie(validReferenceModeSchedule, validStoModeSchedule) = extractValidModeSchedulePair(initTime, finalTime, 
                                                                                            referenceModeSchedule, stoModeSchedule);
  return std::make_pair(std::move(validReferenceModeSchedule.eventTimes), std::move(validStoModeSchedule.eventTimes));
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

std::pair<ModeSchedule, ModeSchedule> extractValidModeSchedulePair(scalar_t initTime, scalar_t finalTime, 
                                                                   const ModeSchedule& referenceModeSchedule,
                                                                   const ModeSchedule& stoModeSchedule) {
  assert(stoModeSchedule.eventTimes.size() == referenceModeSchedule.eventTimes.size());
  assert(stoModeSchedule.modeSequence.size() == referenceModeSchedule.modeSequence.size());
  ModeSchedule validStoModeSchedule, validReferenceModeSchedule;
  size_t lastIndex = 0;
  for (size_t i = 0 ; i < referenceModeSchedule.eventTimes.size(); ++i) {
    lastIndex = i;
    if (referenceModeSchedule.eventTimes[i] > initTime && referenceModeSchedule.eventTimes[i] < finalTime) {
      validReferenceModeSchedule.eventTimes.push_back(referenceModeSchedule.eventTimes[i]);
      validReferenceModeSchedule.modeSequence.push_back(referenceModeSchedule.modeSequence[i]);

      validStoModeSchedule.eventTimes.push_back(stoModeSchedule.eventTimes[i]);
      validStoModeSchedule.modeSequence.push_back(stoModeSchedule.modeSequence[i]);
    }
    else if (referenceModeSchedule.eventTimes[i] >= finalTime) break;
  }
  validReferenceModeSchedule.modeSequence.push_back(referenceModeSchedule.modeSequence[lastIndex]);
  validStoModeSchedule.modeSequence.push_back(stoModeSchedule.modeSequence[lastIndex]);
  return std::make_pair(std::move(validReferenceModeSchedule), std::move(validStoModeSchedule));
}

}  // namespace ocs2

