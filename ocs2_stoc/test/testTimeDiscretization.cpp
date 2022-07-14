#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <gtest/gtest.h>

#include <ocs2_stoc/TimeDiscretization.h>
#include <ocs2_core/NumericTraits.h>

using namespace ocs2;
using namespace stoc;


TEST(testTimeDiscretization, testTimeStepsWithoutSwitches) {
  const scalar_t initTime  = 0.5;
  const scalar_t finalTime = 5.0; 
  const scalar_t dt = 0.03;
  const scalar_array_t eventTimes = {};
  const auto timeDiscretization = multiPhaseTimeDiscretization(initTime, 
                                                               finalTime, dt,
                                                               eventTimes);
  // Verifies initial and final times.
  EXPECT_NEAR(timeDiscretization.front().time, initTime, 
              numeric_traits::weakEpsilon<scalar_t>());
  EXPECT_NEAR(timeDiscretization.back().time, finalTime, 
              numeric_traits::weakEpsilon<scalar_t>());
  // Verifies that the all time intervals are the same.
  const size_t finalStage = timeDiscretization.size() - 1;
  const scalar_t dtPhase = timeDiscretization[1].time - timeDiscretization[0].time;
  for (size_t i=0; i<=finalStage; ++i) {
    EXPECT_NEAR(timeDiscretization[i].time, initTime+i*dtPhase, 
                numeric_traits::weakEpsilon<scalar_t>());
  }
  // Verifies that the size is the same as timeDiscretizationWithEvents.
  const auto timeDiscretizationEvent = timeDiscretizationWithEvents(initTime, 
                                                                    finalTime, dt,
                                                                    eventTimes);
  EXPECT_EQ(timeDiscretizationEvent.size(), timeDiscretization.size());
}

TEST(testTimeDiscretization, testTimeStepsWithSwitches) {
  const scalar_t initTime  = 0.5;
  const scalar_t finalTime = 5.0; 
  const scalar_t dt = 0.06;
  const scalar_array_t eventTimes = {1.0, 1.1, 1.5, 4.0};
  const auto timeDiscretization = multiPhaseTimeDiscretization(initTime, finalTime, dt, eventTimes);
  size_array_t preEventStages;
  for (size_t i=0; i<timeDiscretization.size(); ++i) {
    if (timeDiscretization[i].event == AnnotatedTime::Event::PreEvent) {
      preEventStages.push_back(i);
    }
  }
  // Verifies initial, final, and event times.
  EXPECT_NEAR(timeDiscretization.front().time, initTime, 
              numeric_traits::weakEpsilon<scalar_t>());
  EXPECT_NEAR(timeDiscretization.back().time, finalTime, 
              numeric_traits::weakEpsilon<scalar_t>());
  for (size_t event=0; event<eventTimes.size(); ++event) {
    EXPECT_NEAR(timeDiscretization[preEventStages[event]].time, eventTimes[event],
                numeric_traits::weakEpsilon<scalar_t>());
    EXPECT_NEAR(timeDiscretization[preEventStages[event]+1].time, eventTimes[event],
                numeric_traits::weakEpsilon<scalar_t>());
  }
  // Verifies that the all time intervals are the same.
  scalar_t lastEventTime = initTime; 
  size_t lastPostEventStage = 0;
  for (size_t phase=0; phase<eventTimes.size(); ++phase) {
    const auto dtPhase = timeDiscretization[lastPostEventStage+1].time 
                          - timeDiscretization[lastPostEventStage].time;
    for (size_t stage=lastPostEventStage; stage<=preEventStages[phase]; ++stage) {
      EXPECT_NEAR(lastEventTime+(stage-lastPostEventStage)*dtPhase, 
                  timeDiscretization[stage].time,
                  numeric_traits::weakEpsilon<scalar_t>());
    }
    lastEventTime = eventTimes[phase];
    lastPostEventStage = preEventStages[phase] + 1;
  }
  const auto dtPhase = timeDiscretization[lastPostEventStage+1].time 
                        - timeDiscretization[lastPostEventStage].time;
  const auto finalPhase = eventTimes.size();
  const size_t finalStage = timeDiscretization.size() - 1;
  for (size_t stage=lastPostEventStage; stage<=finalStage; ++stage) {
    EXPECT_NEAR(lastEventTime+(stage-lastPostEventStage)*dtPhase, 
                timeDiscretization[stage].time, 
                numeric_traits::weakEpsilon<scalar_t>());
  }
}

TEST(testTimeDiscretization, testTimeStepsWithSwitches2) {
  const scalar_t initTime  = 1.2;
  const scalar_t finalTime = 5.0; 
  const scalar_t dt = 0.06;
  const scalar_array_t eventTimes = {1.0, 1.1, 1.5, 3.3, 4.0, 5.5, 6.0};
  const auto timeDiscretization = multiPhaseTimeDiscretization(initTime, finalTime, dt, eventTimes);
  size_array_t preEventStages;
  for (size_t i=0; i<timeDiscretization.size(); ++i) {
    if (timeDiscretization[i].event == AnnotatedTime::Event::PreEvent) {
      preEventStages.push_back(i);
    }
  }
  // Verifies initial, final, and event times.
  EXPECT_NEAR(timeDiscretization.front().time, initTime, 
              numeric_traits::weakEpsilon<scalar_t>());
  EXPECT_NEAR(timeDiscretization.back().time, finalTime, 
              numeric_traits::weakEpsilon<scalar_t>());
  // the first two events are skipped as the initTime are greater than these event times
  // the last two events are skipped as the finalTime are less than these event times
  for (size_t event=0; event<eventTimes.size()-4; ++event) {
    EXPECT_NEAR(timeDiscretization[preEventStages[event]].time, eventTimes[event+2],
                numeric_traits::weakEpsilon<scalar_t>());
    EXPECT_NEAR(timeDiscretization[preEventStages[event]+1].time, eventTimes[event+2],
                numeric_traits::weakEpsilon<scalar_t>());
  }
  // Verifies that the all time intervals are the same.
  scalar_t lastEventTime = initTime; 
  size_t lastPostEventStage = 0;
  for (size_t phase=0; phase<eventTimes.size()-4; ++phase) {
    const auto dtPhase = timeDiscretization[lastPostEventStage+1].time 
                          - timeDiscretization[lastPostEventStage].time;
    for (size_t stage=lastPostEventStage; stage<=preEventStages[phase]; ++stage) {
      EXPECT_NEAR(lastEventTime+(stage-lastPostEventStage)*dtPhase, 
                  timeDiscretization[stage].time,
                  numeric_traits::weakEpsilon<scalar_t>());
    }
    lastEventTime = eventTimes[phase+2];
    lastPostEventStage = preEventStages[phase] + 1;
  }
  const auto dtPhase = timeDiscretization[lastPostEventStage+1].time 
                        - timeDiscretization[lastPostEventStage].time;
  const auto finalPhase = eventTimes.size();
  const size_t finalStage = timeDiscretization.size() - 1;
  for (size_t stage=lastPostEventStage; stage<=finalStage; ++stage) {
    EXPECT_NEAR(lastEventTime+(stage-lastPostEventStage)*dtPhase, 
                timeDiscretization[stage].time, 
                numeric_traits::weakEpsilon<scalar_t>());
  }
}