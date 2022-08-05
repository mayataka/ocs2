#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>

namespace ocs2 {

/**
 * Gets the number of valid switching times.
 * @param [in] initTime: The initial time.
 * @param [in] finalTime: The final time.
 * @param [in] modeSchedule: The mode schedule.
 * @return Number of valid switching times.
 */
size_t getNumValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

/**
 * Extracts the valid switching time indices.
 * @param [in] initTime: The initial time.
 * @param [in] finalTime: The final time.
 * @param [in] modeSchedule: The mode schedule.
 * @return Valid switching time indices in the input mode schedule.
 */
size_array_t extractValidSwitchingTimeIndices(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

/**
 * Extracts the valid switching times.
 * @param [in] initTime: The initial time.
 * @param [in] finalTime: The final time.
 * @param [in] modeSchedule: The mode schedule.
 * @return Valid switching times in the input mode schedule.
 */
scalar_array_t extractValidSwitchingTimes(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

/**
 * Extracts the pair of valid switching times.  
 * @param [in] initTime: The initial time.
 * @param [in] finalTime: The final time.
 * @param [in] referenceModeSchedule: The reference mode schedule. This is uesd to evaluate the validness.
 * @param [in] stoModeSchedule: The sto mode schedule.
 * @return Pair of valid switching times. First: valid switching times in reference mode schedule and Second: valid ones in sto mode schedule.
 */
std::pair<scalar_array_t, scalar_array_t> extractValidSwitchingTimesPair(scalar_t initTime, scalar_t finalTime, 
                                                                         const ModeSchedule& referenceModeSchedule,
                                                                         const ModeSchedule& stoModeSchedule);

/**
 * Extracts the valid mode schedule.
 * @param [in] initTime: The initial time.
 * @param [in] finalTime: The final time.
 * @param [in] modeSchedule: The mode schedule.
 * @return Valid mode schedule in the input mode schedule.
 */
ModeSchedule extractValidModeSchedule(scalar_t initTime, scalar_t finalTime, const ModeSchedule& modeSchedule);

/**
 * Extracts the pair of valid mode schedule.  
 * @param [in] initTime: The initial time.
 * @param [in] finalTime: The final time.
 * @param [in] referenceModeSchedule: The reference mode schedule. This is uesd to evaluate the validness.
 * @param [in] stoModeSchedule: The reference mode schedule.
 * @return Pair of valid mode schedules. First: valid mode schedule in reference mode schedule and Second: valid ones in sto mode schedule.
 */
std::pair<ModeSchedule, ModeSchedule> extractValidModeSchedulePair(scalar_t initTime, scalar_t finalTime, 
                                                                   const ModeSchedule& referenceModeSchedule,
                                                                   const ModeSchedule& stoModeSchedule);

}  // namespace ocs2
