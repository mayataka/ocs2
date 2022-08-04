#include <ocs2_stoc/STOC_HelperFunctions.h>

namespace ocs2 {

std::vector<bool> getIsStoEnabled(const ModeSchedule& modeSchedule, const std::unordered_map<size_t, bool>& isStoEnabledInMode) {
  if (modeSchedule.eventTimes.empty()) {
    return std::vector<bool>({}); 
  } else if (modeSchedule.eventTimes.size() == 1) {
    return std::vector<bool>({getisStoEnabledInMode(isStoEnabledInMode, modeSchedule.modeSequence.front())
                              || getisStoEnabledInMode(isStoEnabledInMode, modeSchedule.modeSequence.back())}); 
  } else {
    std::vector<bool> isStoEnabled;
    isStoEnabled.push_back(getisStoEnabledInMode(isStoEnabledInMode, modeSchedule.modeSequence.front()));
    for (size_t phase=1; phase<modeSchedule.eventTimes.size()-1; ++phase) {
      isStoEnabled.push_back(getisStoEnabledInMode(isStoEnabledInMode, modeSchedule.modeSequence[phase-1])
                            || getisStoEnabledInMode(isStoEnabledInMode, modeSchedule.modeSequence[phase]));
    }
    isStoEnabled.push_back(getisStoEnabledInMode(isStoEnabledInMode, modeSchedule.modeSequence.back()));
    return isStoEnabled;
  }
}

bool getisStoEnabledInMode(const std::unordered_map<size_t, bool>& isStoEnabledInMode, size_t mode) {
  if (isStoEnabledInMode.find(mode) == isStoEnabledInMode.end()) return false;
  else return isStoEnabledInMode.at(mode);
}

}  // namespace ocs2
