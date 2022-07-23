#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>

namespace ocs2 {

scalar_array_t extractValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

ModeSchedule extractValidModeSchedule(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

}  // namespace ocs2
