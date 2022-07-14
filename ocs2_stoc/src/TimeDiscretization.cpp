#include <ocs2_stoc/TimeDiscretization.h>

#include <cassert>
#include <string>

namespace ocs2 {
namespace stoc {

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

}  // namespace stoc 
}  // namespace ocs2


namespace ocs2 {

std::string EventName(const AnnotatedTime::Event& event) {
  switch (event) {
    case AnnotatedTime::Event::None:
      return "None";
      break;
    case AnnotatedTime::Event::PreEvent:
      return "PreEvent";
      break;
    case AnnotatedTime::Event::PostEvent:
      return "PostEvent";
      break;
    default:
      return "";
      break;
    }
}

std::ostream& operator<<(std::ostream& stream, const std::vector<AnnotatedTime>& timeDiscretization) {
  for (const auto& e : timeDiscretization) {
    stream << "time: " << e.time << ",  event type: " << EventName(e.event) << "\n";
  }
  return stream;
}

}  // namespace ocs2