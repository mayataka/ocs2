#include "ocs2_stoc/reference/STOC_ModeSchedule.h"

#include <ocs2_core/misc/Display.h>
#include <ocs2_core/misc/Lookup.h>
#include <ocs2_core/misc/Numerics.h>

namespace ocs2 {
namespace stoc {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
STOC_ModeSchedule::STOC_ModeSchedule(std::vector<scalar_t> eventTimesInput, std::vector<size_t> modeSequenceInput, std::vector<bool> isStoEnabledInput)
    : eventTimes(std::move(eventTimesInput)), modeSequence(std::move(modeSequenceInput)), isStoEnabled(std::move(isStoEnabledInput)) {
  assert(!modeSequence.empty());
  assert(eventTimes.size() + 1 == modeSequence.size());
  assert(eventTimes.size() == isStoEnabled.size());
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
size_t STOC_ModeSchedule::modeAtTime(scalar_t time) const {
  const auto ind = lookup::findIndexInTimeArray(eventTimes, time);
  return modeSequence[ind];
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void swap(STOC_ModeSchedule& lh, STOC_ModeSchedule& rh) {
  lh.eventTimes.swap(rh.eventTimes);
  lh.modeSequence.swap(rh.modeSequence);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::ostream& operator<<(std::ostream& stream, const STOC_ModeSchedule& modeSchedule) {
  stream << "event times:   {" << toDelimitedString(modeSchedule.eventTimes) << "}\n";
  stream << "mode sequence: {" << toDelimitedString(modeSchedule.modeSequence) << "}\n";
  stream << "is STO enabled: {" << toDelimitedString(modeSchedule.isStoEnabled) << "}\n";
  return stream;
}

}  // namespace stoc
}  // namespace ocs2
