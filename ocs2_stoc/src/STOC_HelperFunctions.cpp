#include <ocs2_stoc/STOC_HelperFunctions.h>

namespace ocs2 {

std::vector<bool> getIsStoEnabled(const ModeSchedule& modeSchedule, const std::unordered_map<size_t, bool>& isStoEnabledMode) {
  if (modeSchedule.eventTimes.empty()) {
    return std::vector<bool>({}); 
  } else if (modeSchedule.eventTimes.size() == 1) {
    return std::vector<bool>({getIsStoEnabledMode(isStoEnabledMode, modeSchedule.modeSequence.front())
                              || getIsStoEnabledMode(isStoEnabledMode, modeSchedule.modeSequence.back())}); 
  } else {
    std::vector<bool> isStoEnabled;
    isStoEnabled.push_back(getIsStoEnabledMode(isStoEnabledMode, modeSchedule.modeSequence.front()));
    for (size_t phase=1; phase<modeSchedule.eventTimes.size()-1; ++phase) {
      isStoEnabled.push_back(getIsStoEnabledMode(isStoEnabledMode, modeSchedule.modeSequence[phase-1])
                            || getIsStoEnabledMode(isStoEnabledMode, modeSchedule.modeSequence[phase]));
    }
    isStoEnabled.push_back(getIsStoEnabledMode(isStoEnabledMode, modeSchedule.modeSequence.back()));
    return isStoEnabled;
  }
}

bool getIsStoEnabledMode(const std::unordered_map<size_t, bool>& isStoEnabledMode, size_t mode) {
  if (isStoEnabledMode.find(mode) == isStoEnabledMode.end()) return false;
  else return isStoEnabledMode.at(mode);
}

}  // namespace ocs2
