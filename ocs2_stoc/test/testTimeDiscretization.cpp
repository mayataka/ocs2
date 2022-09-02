#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <gtest/gtest.h>

#include <ocs2_stoc/TimeDiscretization.h>
#include <ocs2_core/NumericTraits.h>
#include <ocs2_core/reference/ModeSchedule.h>

using namespace ocs2;


TEST(testTimeDiscretization, testTimeStepsWithoutSwitches) {
  const scalar_t initTime  = 0.5;
  const scalar_t finalTime = 5.0; 
  const scalar_t dt = 0.03;
  const scalar_array_t eventTimes = {};
  const auto timeDiscretization = multiPhaseTimeDiscretization(initTime, finalTime, dt, eventTimes);
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
  const auto timeDiscretizationEvent = timeDiscretizationWithEvents(initTime, finalTime, dt, eventTimes);
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

TEST(testTimeDiscretization, testCastEvent) {
  EXPECT_EQ(castEvent(AnnotatedTime::Event::None), Grid::Event::None);
  EXPECT_EQ(castEvent(AnnotatedTime::Event::PreEvent), Grid::Event::PreEvent);
  EXPECT_EQ(castEvent(AnnotatedTime::Event::PostEvent), Grid::Event::PostEvent);
}

TEST(testTimeDiscretization, testMultiPhaseTimeDiscretizationGrid) {
  const scalar_t initTime  = 0.5;
  const scalar_t finalTime = 5.0; 
  const scalar_t dt = 0.06;
  const scalar_t min_dt = 10.0 * numeric_traits::limitEpsilon<scalar_t>();
  const scalar_array_t eventTimesInput = {1.0, 1.1, 1.5, 4.0};
  const size_array_t modeScheduleInput = {1, 0, 1, 2, 0};
  const auto modeSchedule = ModeSchedule(eventTimesInput, modeScheduleInput);
  const auto grids = multiPhaseTimeDiscretizationGrid(initTime, finalTime, dt, modeSchedule);
  const auto timeDiscretization = multiPhaseTimeDiscretization(initTime, finalTime, dt, modeSchedule.eventTimes);
  size_t phase = 0;
  for (size_t i=0; i<grids.size(); ++i) {
    EXPECT_EQ(grids[i].time, timeDiscretization[i].time);
    EXPECT_EQ(grids[i].event, castEvent(timeDiscretization[i].event));
    size_t mode = modeSchedule.modeAtTime(grids[i].time);
    if (grids[i].event == Grid::Event::PostEvent) {
      mode = modeSchedule.modeAtTime(grids[i].time + min_dt);
      ++phase;
    }
    EXPECT_EQ(grids[i].mode, mode);
    EXPECT_EQ(grids[i].phase, phase);
  }
  std::cout << grids << std::endl;

  const auto numGrids = getNumGrids(grids);
  size_array_t numGridsRef(numGrids.size(), 0);
  for (size_t i=0; i<grids.size()-1; ++i) {
    if (grids[i].event != Grid::Event::PreEvent) {
      ++numGridsRef[grids[i].phase];
    }
  }
  for (size_t i=0; i<numGrids.size(); ++i) {
    EXPECT_EQ(numGrids[i], numGridsRef[i]);
  }

  for (int i=0; i<numGrids.size(); ++i) {
    std::cout << "numGrids[i]: " <<  numGrids[i] << std::endl;
  }
}

TEST(testTimeDiscretization, testMultiPhaseTimeDiscretizationGridWithSTO) {
  const scalar_t initTime  = 1.1;
  const scalar_t finalTime = 5.0; 
  const scalar_t dt = 0.06;
  const scalar_t min_dt = 10.0 * numeric_traits::limitEpsilon<scalar_t>();
  const scalar_array_t eventTimesInput = {1.0, 1.1, 1.5, 4.0, 5.5};
  const size_array_t modeScheduleInput = {1, 0, 1, 2, 0, 1};
  const auto modeSchedule = ModeSchedule(eventTimesInput, modeScheduleInput);
  const std::unordered_map<size_t, bool> isStoEnabledInMode = {{0, true}, {1, false}, {2, true}, };
  const auto grids = multiPhaseTimeDiscretizationGrid(initTime, finalTime, dt, modeSchedule, isStoEnabledInMode);
  const auto timeDiscretization = multiPhaseTimeDiscretization(initTime, finalTime, dt, modeSchedule.eventTimes);
  size_t phase = 0;
  for (size_t i=0; i<grids.size(); ++i) {
    EXPECT_EQ(grids[i].time, timeDiscretization[i].time);
    if (i > 0) {
      EXPECT_EQ(grids[i].event, castEvent(timeDiscretization[i].event)); // event whose evenTime is initTime is skipped
    }
    size_t mode = modeSchedule.modeAtTime(grids[i].time);
    if (grids[i].event == Grid::Event::PostEvent && i > 0) { // event whose evenTime is initTime is skipped
      mode = modeSchedule.modeAtTime(grids[i].time + min_dt);
      ++phase;
    }
    EXPECT_EQ(grids[i].mode, mode);
    EXPECT_EQ(grids[i].phase, phase);
  }
  std::cout << grids << std::endl;
}

TEST(testTimeDiscretization, testUpdateTimeIntervals) {
  const scalar_t initTime  = 0.5;
  const scalar_t finalTime = 5.0; 
  const scalar_t dt = 0.06;
  const scalar_t min_dt = 10.0 * numeric_traits::limitEpsilon<scalar_t>();
  const scalar_array_t eventTimesInput = {1.0, 1.1, 1.5, 4.0};
  const size_array_t modeScheduleInput = {1, 0, 1, 2, 0};
  const auto modeSchedule = ModeSchedule(eventTimesInput, modeScheduleInput);
  const auto gridsRef = multiPhaseTimeDiscretizationGrid(initTime, finalTime, dt, modeSchedule);
  const scalar_array_t switchingTimeDirections = {0.2, -0.3, 1.0, 0.3};
  auto grids = gridsRef;
  updateTimeIntervals(initTime, finalTime, switchingTimeDirections, grids);
  for (size_t i=0; i<grids.size(); ++i) {
    EXPECT_EQ(grids[i].event, gridsRef[i].event);
    EXPECT_EQ(grids[i].mode, gridsRef[i].mode);
    EXPECT_EQ(grids[i].phase, gridsRef[i].phase);
    EXPECT_EQ(grids[i].sto, gridsRef[i].sto);
    EXPECT_EQ(grids[i].stoNext, gridsRef[i].stoNext);
    EXPECT_EQ(grids[i].stoNextNext, gridsRef[i].stoNextNext);
  }
  const auto numGrids = getNumGrids(gridsRef);
  scalar_array_t updatedEventTimes = modeSchedule.eventTimes;
  for (size_t i = 0; i < updatedEventTimes.size(); ++i) {
    updatedEventTimes[i] += switchingTimeDirections[i];
  }
  scalar_array_t dtPhase(updatedEventTimes.size() + 1); 
  dtPhase[0] = (updatedEventTimes[0] - initTime) / numGrids[0];
  for (size_t i=1; i<dtPhase.size() - 1; ++i) {
    dtPhase[i] = (updatedEventTimes[i] - updatedEventTimes[i - 1]) / numGrids[i];
  }
  dtPhase.back() = (finalTime - updatedEventTimes.back()) / numGrids.back();
  scalar_t time = initTime;
  for (size_t i=0; i<grids.size(); ++i) {
    EXPECT_NEAR(grids[i].time, time, numeric_traits::weakEpsilon<scalar_t>());
    if (grids[i].event != Grid::Event::PreEvent) {
      time += dtPhase[grids[i].phase];
    }
  }
}