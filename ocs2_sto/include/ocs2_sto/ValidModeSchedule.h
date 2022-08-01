#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>

namespace ocs2 {

size_t getNumValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

size_array_t extractValidSwitchingTimeIndices(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

scalar_array_t extractValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

std::pair<scalar_array_t, scalar_array_t> extractValidSwitchingTimesPair(scalar_t initTime, scalar_t finalTime, 
                                                                         const ModeSchedule& referenceModeSchedule,
                                                                         const ModeSchedule& stoModeSchedule);

ModeSchedule extractValidModeSchedule(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

std::pair<ModeSchedule, ModeSchedule> extractValidModeSchedulePair(scalar_t initTime, scalar_t finalTime, 
                                                                   const ModeSchedule& referenceModeSchedule,
                                                                   const ModeSchedule& stoModeSchedule);

}  // namespace ocs2
