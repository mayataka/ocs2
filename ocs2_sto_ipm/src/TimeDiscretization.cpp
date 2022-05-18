#include <ocs2_sto_ipm/TimeDiscretization.h>

#include <cassert>

namespace ocs2 {
namespace sto_ipm {

std::vector<AnnotatedTime> multiPhaseTimeDiscretization(scalar_t initTime, scalar_t finalTime, scalar_t dt,
                                                        const scalar_array_t& eventTimes, scalar_t dt_min) {
  assert(dt > 0);
  assert(finalTime > initTime);
  std::vector<AnnotatedTime> timeDiscretization
      = timeDiscretizationWithEvents(initTime, finalTime, dt, eventTimes, dt_min);
  const auto numPhases = eventTimes.size() + 1;
  scalar_t lastEventTime = initTime;
  size_t lastPostEventGrid = 0;
  for (size_t phase=0; phase<numPhases-1; ++phase) {
    size_t nextPreEventGrid = 0;
    for (size_t j=lastPostEventGrid; j<timeDiscretization.size(); ++j) {
      if (timeDiscretization[j].event == AnnotatedTime::Event::PreEvent) {
        nextPreEventGrid = j;
        break;
      }
    }
    const auto numGridInPhase = nextPreEventGrid - lastPostEventGrid;
    const auto dtPhase = (eventTimes[phase] - lastEventTime) / numGridInPhase;
    for (size_t j=lastPostEventGrid; j<nextPreEventGrid; ++j) {
      timeDiscretization[j].time = lastEventTime + (j-lastPostEventGrid) * dtPhase;
    }
    lastPostEventGrid = nextPreEventGrid + 1;
    lastEventTime = timeDiscretization[lastPostEventGrid].time;
  }
  const size_t finalPhase = numPhases - 1;
  const size_t finalGrid = timeDiscretization.size() - 1;
  const auto numGridInPhase = finalGrid - lastPostEventGrid;
  const auto dtPhase = (finalTime - lastEventTime) / numGridInPhase;
  for (size_t j=lastPostEventGrid; j<finalGrid; ++j) {
    timeDiscretization[j].time = lastEventTime + (j-lastPostEventGrid) * dtPhase;
  }
  return timeDiscretization;
}

}  // namespace sto_ipm 
}  // namespace ocs2
