#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>

namespace ocs2 {

scalar_array_t extractValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

std::pair<scalar_array_t, scalar_array_t> extractValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, 
                                                                     const ModeSchedule& stoModeSchedule,
                                                                     const ModeSchedule& referenceModeSchedule);

ModeSchedule extractValidModeSchedule(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

std::pair<ModeSchedule, ModeSchedule> extractValidModeSchedule(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule,
                                                               const ModeSchedule& referenceModeSchedule);

}  // namespace ocs2
