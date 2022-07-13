#include <ocs2_stoc/DiscreteTimeModeSchedule.h>

#include <ocs2_core/misc/Display.h>
#include <ocs2_core/misc/Lookup.h>
#include <ocs2_core/misc/Numerics.h>

#include <cassert>
#include <iostream>

namespace ocs2 {
namespace stoc {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
DiscreteTimeModeSchedule::DiscreteTimeModeSchedule(std::vector<size_t> modeSequenceInput, std::vector<bool> isStoEnabledInput)
    : modeSequence(std::move(modeSequenceInput)), phaseSequence({}), isStoEnabled() {
  assert(!modeSequence.empty());
  auto modePrev = modeSequence.front();
  size_t phase = 0;
  for (const auto mode : modeSequence) {
    if (mode != modePrev) {
      ++phase;
    }
    phaseSequence.push_back(phase);
    isStoEnabled.push_back(isStoEnabledInput[phase]);
    modePrev = mode;
  }
  if (!isStoEnabledInput.empty()) {
    assert(isStoEnabledInput.size() == phaseSequence.back());
    for (const auto phase : phaseSequence) {
      isStoEnabled.push_back(isStoEnabledInput[phase]);
    }
  }
  else {
    isStoEnabled = std::vector<bool>(modeSequence.size()-1, false);
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
DiscreteTimeModeSchedule::DiscreteTimeModeSchedule(const std::vector<AnnotatedTime>& timeDiscretization,
                                                   const STOC_ModeSchedule& modeScheduleReference) 
    : modeSequence({}), phaseSequence({}), isStoEnabled({}) {
  assert(!modeSequence.empty());
  auto modePrev = modeScheduleReference.modeAtTime(timeDiscretization.front().time);
  size_t phase = 0;
  for (const auto& e : timeDiscretization) {
    const auto mode = modeScheduleReference.modeAtTime(e.time);
    if (mode != modePrev) {
      ++phase;
    }
    modeSequence.push_back(mode);
    phaseSequence.push_back(phase);
    isStoEnabled.push_back(modeScheduleReference.isStoEnabled[phase]);
    modePrev = mode;
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
size_t DiscreteTimeModeSchedule::modeAtTimeStage(size_t timeStage) const {
  assert(timeStage >= 0);
  assert(timeStage < modeSequence.size());
  return modeSequence[timeStage];
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
size_t DiscreteTimeModeSchedule::phaseAtTimeStage(size_t timeStage) const {
  assert(timeStage >= 0);
  assert(timeStage < phaseSequence.size());
  return phaseSequence[timeStage];
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
bool DiscreteTimeModeSchedule::isStoEnabledAtPhase(size_t phase) const {
  if (phase >= phaseSequence.size()) {
    return false;
  }
  else {
    return isStoEnabled[phase];
  }
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void swap(DiscreteTimeModeSchedule& lh, DiscreteTimeModeSchedule& rh) {
  lh.modeSequence.swap(rh.modeSequence);
  lh.phaseSequence.swap(rh.phaseSequence);
  lh.isStoEnabled.swap(rh.isStoEnabled);
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::ostream& operator<<(std::ostream& stream, const DiscreteTimeModeSchedule& modeSchedule) {
  stream << "mode sequence: {" << toDelimitedString(modeSchedule.modeSequence) << "}\n";
  stream << "phase sequence: {" << toDelimitedString(modeSchedule.phaseSequence) << "}\n";
  stream << "is STO enabled at each phase:   {" << toDelimitedString(modeSchedule.isStoEnabled) << "}\n";
  return stream;
}

}  // namespace stoc
}  // namespace ocs2
