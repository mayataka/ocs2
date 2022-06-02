#include <iostream>
#include <typeinfo>

#include <gtest/gtest.h>

#include <ocs2_stoc/reference/DiscreteTimeModeSchedule.h>

using namespace ocs2;
using namespace stoc;

TEST(testDiscreteTimeModeSchedule, testDiscreteTimeModeSchedule) {
  // 5 phases, 4 switches
  const size_array_t modeSequenceRef = {1, 1, 1, 0, 0, 0, 0, 2, 2, 1, 0, 0, 0}; 
  const size_array_t modeSequenceInput = modeSequenceRef;
  const std::vector<bool> isStoEnabledRef = {true, false, false, true};
  const std::vector<bool> isStoEnabledInput = isStoEnabledRef;

  DiscreteTimeModeSchedule modeSchedule(modeSequenceInput, isStoEnabledInput);

  for (int i=0; i<modeSequenceRef.size(); ++i) {
    EXPECT_EQ(modeSchedule.modeAtTimeStage(i), modeSequenceRef[i]);
  }
  for (int i=0; i<isStoEnabledRef.size(); ++i) {
    EXPECT_EQ(modeSchedule.isStoEnabledAtPhase(i), isStoEnabledInput[i]);
  }
  // mode 1: 3 time stages 
  for (int i=0; i<3; ++i) {
    EXPECT_EQ(modeSchedule.modeAtTimeStage(i), 1);
    EXPECT_EQ(modeSchedule.phaseAtTimeStage(i), 0);
  }
  // mode 0: 4 times stages
  for (int i=3; i<7; ++i) {
    EXPECT_EQ(modeSchedule.modeAtTimeStage(i), 0);
    EXPECT_EQ(modeSchedule.phaseAtTimeStage(i), 1);
  }
  // mode 2: 2 times stages
  for (int i=7; i<9; ++i) {
    EXPECT_EQ(modeSchedule.modeAtTimeStage(i), 2);
    EXPECT_EQ(modeSchedule.phaseAtTimeStage(i), 2);
  }
  // mode 1: 1 times stages
  EXPECT_EQ(modeSchedule.modeAtTimeStage(9), 1);
  EXPECT_EQ(modeSchedule.phaseAtTimeStage(9), 3);
  // mode 0: 3 times stages
  for (int i=10; i<13; ++i) {
    EXPECT_EQ(modeSchedule.modeAtTimeStage(i), 0);
    EXPECT_EQ(modeSchedule.phaseAtTimeStage(i), 4);
  }
}
