#include <ocs2_stoc/reference/DiscreteTimeModeSchedule.h>

#include <ocs2_core/misc/Display.h>
#include <ocs2_core/misc/Lookup.h>
#include <ocs2_core/misc/Numerics.h>

#include <cassert>

namespace ocs2 {
namespace stoc {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
DiscreteTimeModeSchedule::DiscreteTimeModeSchedule(std::vector<size_t> modeSequenceInput, 
                                                   std::vector<bool> isStoEnabled)
    : modeSequence(std::move(modeSequenceInput)), isStoEnabled(std::move(isStoEnabled)) {
  assert(!modeSequence.empty());
  phaseSequence.clear();
  size_t modePrev = modeSequence[0];
  size_t phase = 0;
  for (const auto mode : modeSequence) {
    if (mode != modePrev) {
      ++phase;
      modePrev = mode;
    }
    phaseSequence.push_back(phase);
  }
  assert(phase == isStoEnabled.size());
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
