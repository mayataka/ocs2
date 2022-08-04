#pragma once

#include <unordered_map>

#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>

namespace ocs2 {

std::vector<bool> getIsStoEnabled(const ModeSchedule& modeSchedule, const std::unordered_map<size_t, bool>& isStoEnabledInMode);

bool getisStoEnabledInMode(const std::unordered_map<size_t, bool>& isStoEnabledInMode, size_t mode);

}  // namespace ocs2
