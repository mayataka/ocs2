#include "ocs2_stoc/TimeDiscretization.h"

#include <cassert>
#include <string>

namespace ocs2 {

scalar_t getInterpolationTime(const Grid& grid) {
  return grid.time + numeric_traits::limitEpsilon<scalar_t>();
}

scalar_t getIntervalStart(const Grid& start) {
  scalar_t adaptedStart = start.time;
  if (start.event == Grid::Event::PostEvent) {
    adaptedStart += numeric_traits::weakEpsilon<scalar_t>();
  }
  return adaptedStart;
}

scalar_t getIntervalEnd(const Grid& end) {
  scalar_t adaptedEnd = end.time;
  if (end.event == Grid::Event::PreEvent) {
    adaptedEnd -= numeric_traits::weakEpsilon<scalar_t>();
  }
  return adaptedEnd;
}

scalar_t getIntervalDuration(const Grid& start, const Grid& end) {
  return getIntervalEnd(end) - getIntervalStart(start);
}

std::vector<Grid> multiPhaseTimeDiscretizationGrid(scalar_t initTime, scalar_t finalTime, scalar_t dt, const ModeSchedule& modeSchedule, 
                                                   const std::unordered_map<size_t, bool>& isStoEnabledInMode, scalar_t dt_min) {
  auto timeDiscretization = multiPhaseTimeDiscretization(initTime, finalTime, dt, modeSchedule.eventTimes, dt_min);
  timeDiscretization.front().event = AnnotatedTime::Event::None; // disable post-event at the initial stage.
  std::vector<Grid> timeDiscretizationGrid;
  timeDiscretizationGrid.reserve(timeDiscretization.size());
  size_t phase = 0;
  for (const auto& e : timeDiscretization) {
    auto mode = modeSchedule.modeAtTime(e.time);
    if (e.event == AnnotatedTime::Event::PostEvent) {
      mode = modeSchedule.modeAtTime(e.time+dt_min); // TODO: test this carefully 
      ++phase;
    } 
    constexpr bool sto = false;
    constexpr bool stoNext = false;
    timeDiscretizationGrid.emplace_back(e.time, mode, phase, castEvent(e.event), sto, stoNext);
  }

  auto extractValidModeSequence = [](const std::vector<Grid>& timeDiscretizationGrid) {
    std::vector<size_t> validModeSequence;
    validModeSequence.reserve(timeDiscretizationGrid.back().phase+1); 
    auto prevPhase = timeDiscretizationGrid.front().phase;
    validModeSequence.push_back(timeDiscretizationGrid.front().mode);
    for (const auto& e : timeDiscretizationGrid) {
      if (e.phase != prevPhase) {
        validModeSequence.push_back(e.mode);
        prevPhase = e.phase;
      }
    }
    return validModeSequence;
  };
  const std::vector<size_t> validModeSequence = extractValidModeSequence(timeDiscretizationGrid);

  auto checkIsStoEnabledInMode = [](const std::unordered_map<size_t, bool>& isStoEnabledInMode, size_t mode) {
    if (isStoEnabledInMode.empty()) return false;
    if (isStoEnabledInMode.find(mode) == isStoEnabledInMode.end()) return false;
    return isStoEnabledInMode.at(mode); 
  };
  std::vector<bool> isStoEnabledInPhase;
  isStoEnabledInPhase.reserve(validModeSequence.size()+2); // maximum size
  for (const auto mode : validModeSequence) {
    isStoEnabledInPhase.push_back(checkIsStoEnabledInMode(isStoEnabledInMode, mode));
  }
  isStoEnabledInPhase.push_back(false);
  isStoEnabledInPhase.push_back(false);
  for (auto& e : timeDiscretizationGrid) {
    e.sto = isStoEnabledInPhase[e.phase];
    if (e.phase+1 < isStoEnabledInPhase.size()) {
      e.stoNext = isStoEnabledInPhase[e.phase+1];
    }
  }
  return timeDiscretizationGrid; 
}

void updateTimeIntervals(scalar_t initTime, scalar_t finalTime, const scalar_array_t& switchingTimeDirections, 
                         std::vector<Grid>& timeDiscretization) {
  if (switchingTimeDirections.empty()) return;

  const int N = static_cast<int>(timeDiscretization.size()) - 1;

  size_t validEventIndex = 0;
  for (int i=0; i<N; ++i) {
    if (timeDiscretization[i].event == Grid::Event::PreEvent) {
      timeDiscretization[i].time += switchingTimeDirections[validEventIndex];
    } else if (timeDiscretization[i].event == Grid::Event::PostEvent) {
      timeDiscretization[i].time += switchingTimeDirections[validEventIndex];
      ++validEventIndex;
    }
  }

  int lastPostEventGrid = 0;
  scalar_t lastEventTime = initTime;
  for (int i=0; i<N; ++i) {
    if (timeDiscretization[i].event == Grid::Event::PreEvent) {
      const int numGrids = i - lastPostEventGrid;
      const scalar_t dt = (timeDiscretization[i].time - lastEventTime) / numGrids;
      for (int j=lastPostEventGrid+1; j<i; ++j) {
        timeDiscretization[j].time = lastEventTime + (j - lastPostEventGrid) * dt;
      }
    } else if (timeDiscretization[i].event == Grid::Event::PostEvent) {
      lastPostEventGrid = i;
      lastEventTime = timeDiscretization[i].time;
    }
  }
  const int numGrids = N - lastPostEventGrid;
  const scalar_t dt = (timeDiscretization[N].time - lastEventTime) / numGrids;
  for (int j=lastPostEventGrid+1; j<N; ++j) {
    timeDiscretization[j].time = lastEventTime + (j - lastPostEventGrid) * dt;
  }
  timeDiscretization[N].time = finalTime;
}

std::vector<AnnotatedTime> multiPhaseTimeDiscretization(scalar_t initTime, scalar_t finalTime, scalar_t dt,
                                                        const scalar_array_t& eventTimes, scalar_t dt_min) {
  assert(dt > 0);
  assert(finalTime > initTime);
  std::vector<AnnotatedTime> timeDiscretization
      = timeDiscretizationWithEvents(initTime, finalTime, dt, eventTimes, dt_min);
  if (eventTimes.empty()) {
    return timeDiscretization;
  }
  const auto numPhases = eventTimes.size() + 1;
  scalar_t lastEventTime = initTime;
  int lastPostEventGrid = 0;
  int initPhase, finalPhase;
  for (initPhase=0; initPhase<numPhases; ++initPhase) {
    if (eventTimes[initPhase] > initTime) break;
  }
  for (finalPhase=initPhase; finalPhase<numPhases-1; ++finalPhase) {
    if (eventTimes[finalPhase+1] >= finalTime) break;
  }
  if (finalPhase == numPhases-1) --finalPhase;
  for (int phase=initPhase; phase<=finalPhase; ++phase) {
    int nextPreEventGrid = 0;
    for (int j=lastPostEventGrid; j<timeDiscretization.size(); ++j) {
      if (timeDiscretization[j].event == AnnotatedTime::Event::PreEvent) {
        nextPreEventGrid = j;
        break;
      }
    }
    const auto numGridInPhase = nextPreEventGrid - lastPostEventGrid;
    const auto dtPhase = (eventTimes[phase] - lastEventTime) / numGridInPhase;
    for (int j=lastPostEventGrid; j<nextPreEventGrid; ++j) {
      timeDiscretization[j].time = lastEventTime + (j-lastPostEventGrid) * dtPhase;
    }
    lastPostEventGrid = nextPreEventGrid + 1;
    lastEventTime = timeDiscretization[lastPostEventGrid].time;
  }
  const int finalGrid = timeDiscretization.size() - 1;
  const auto numGridInPhase = finalGrid - lastPostEventGrid;
  const auto dtPhase = (finalTime - lastEventTime) / numGridInPhase;
  for (int j=lastPostEventGrid; j<finalGrid; ++j) {
    timeDiscretization[j].time = lastEventTime + (j-lastPostEventGrid) * dtPhase;
  }
  return timeDiscretization;
}

scalar_t getMaxTimeInterval(const std::vector<Grid>& timeDiscretizationGrid) {
  assert(!timeDiscretizationGrid.empty());
  scalar_t maxTimeInterval = timeDiscretizationGrid[1].time - timeDiscretizationGrid[0].time;
  for (size_t i = 0; i < timeDiscretizationGrid.size() - 1; ++i) {
    if (timeDiscretizationGrid[i].event == Grid::Event::PostEvent) {
      maxTimeInterval = std::max(maxTimeInterval, (timeDiscretizationGrid[i + 1].time - timeDiscretizationGrid[i].time));
    } 
  }
  return maxTimeInterval;
}

size_array_t getNumGrids(const std::vector<Grid>& timeDiscretizationGrid) {
  size_array_t numGrids(timeDiscretizationGrid.back().phase+1, 0);
  for (size_t i=0; i < timeDiscretizationGrid.size()-1; ++i) {
    if (timeDiscretizationGrid[i].event != Grid::Event::PreEvent) {
      ++numGrids[timeDiscretizationGrid[i].phase];
    } 
  }
  return numGrids;
}

std::string toString(const Grid::Event& event) {
  switch (event) {
    case Grid::Event::None:
      return "None";
      break;
    case Grid::Event::PreEvent:
      return "PreEvent";
      break;
    case Grid::Event::PostEvent:
      return "PostEvent";
      break;
    default:
      return "";
      break;
    }
}

std::ostream& operator<<(std::ostream& stream, const std::vector<AnnotatedTime>& timeDiscretization) {
  for (const auto& e : timeDiscretization) {
    stream << "time: " << e.time << ",  event type: " << toString(castEvent(e.event)) << "\n";
  }
  return stream;
}

std::ostream& operator<<(std::ostream& stream, const std::vector<Grid>& timeDiscretizationGrid) {
  for (const auto& e : timeDiscretizationGrid) {
    stream << "time: " << e.time 
           << ", mode: " << e.mode 
           << ", phase: " << e.phase 
           << ", event: " << toString(e.event) 
           << ", sto: " << std::boolalpha << e.sto
           << ", stoNext: " << std::boolalpha << e.stoNext
           << "\n";
  }
  return stream;
}

}  // namespace ocs2