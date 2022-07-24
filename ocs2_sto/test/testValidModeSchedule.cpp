#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <gtest/gtest.h>

#include <ocs2_core/NumericTraits.h>
#include <ocs2_core/reference/ModeSchedule.h>
#include <ocs2_sto/ValidModeSchedule.h>

using namespace ocs2;


TEST(testTimeDiscretization, testTimeDiscretization) {
  const scalar_t initTime  = 0.5;
  const scalar_t finalTime = 5.0; 
  const scalar_array_t eventTimesInput = {1.0, 1.1, 1.5, 4.0, 5.0, 6.0};
  const size_array_t modeScheduleInput = {1, 0, 1, 2, 0, 1, 2};
  const auto modeSchedule = ModeSchedule(eventTimesInput, modeScheduleInput);
  const auto validModeSchedule = extractValidModeSchedule(initTime, finalTime, modeSchedule);
  for (const auto e : validModeSchedule.eventTimes) {
    EXPECT_TRUE(e > initTime);
    EXPECT_TRUE(e < finalTime);
  }
  const auto validSwitchingTimes = extractValidSwitchingTimes(initTime, finalTime, modeSchedule);
  for (int i=0; i<validSwitchingTimes.size(); ++i) {
    EXPECT_DOUBLE_EQ(validSwitchingTimes[i], validModeSchedule.eventTimes[i]);
  }

  const auto referenceModeSchedule = modeSchedule;
  auto stoModeSchedule = modeSchedule;
  for (auto& e : stoModeSchedule.eventTimes) {
    e += 5.0; 
  }

  const auto validModeSchedulePair = extractValidModeSchedule(initTime, finalTime, stoModeSchedule, referenceModeSchedule);
  const auto validStoModeSchedule = validModeSchedulePair.first;
  const auto validReferenceModeSchedule = validModeSchedulePair.second;
  for (const auto e : validReferenceModeSchedule.eventTimes) {
    EXPECT_TRUE(e > initTime);
    EXPECT_TRUE(e < finalTime);
  }
  for (int i=0; i<validStoModeSchedule.eventTimes.size(); ++i) {
    EXPECT_DOUBLE_EQ(validStoModeSchedule.eventTimes[i]-5.0, validReferenceModeSchedule.eventTimes[i]);
  }
  for (int i=0; i<validStoModeSchedule.modeSequence.size(); ++i) {
    EXPECT_EQ(validStoModeSchedule.modeSequence[i], validReferenceModeSchedule.modeSequence[i]);
  }

  const auto validSwitchingTimesPair = extractValidSwitchingTimes(initTime, finalTime, stoModeSchedule, referenceModeSchedule);
  const auto validStoSwitchingTimes = validSwitchingTimesPair.first;
  const auto validReferenceSwitchingTimes = validSwitchingTimesPair.second;
  for (int i=0; i<validStoSwitchingTimes.size(); ++i) {
    EXPECT_DOUBLE_EQ(validStoSwitchingTimes[i], validStoModeSchedule.eventTimes[i]);
  }
  for (int i=0; i<validReferenceSwitchingTimes.size(); ++i) {
    EXPECT_DOUBLE_EQ(validReferenceSwitchingTimes[i], validReferenceModeSchedule.eventTimes[i]);
  }
}