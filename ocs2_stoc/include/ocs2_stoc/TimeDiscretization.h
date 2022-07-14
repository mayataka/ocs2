#pragma once

#include <ocs2_core/NumericTraits.h>
#include <ocs2_core/Types.h>
#include <ocs2_sqp/TimeDiscretization.h>
#include <ocs2_core/reference/ModeSchedule.h>

#include <iostream>

namespace ocs2 {

struct Grid {
  /// Enum defines the type of event that occurs at this time
  enum class Event { None, PreEvent, PostEvent };

  scalar_t time;
  size_t mode;
  size_t phase;
  Event event;
  bool sto;
  bool stoNext;
  bool stoNextNext;

  /** Constructor defaulting to 'None' event  */
  explicit Grid(scalar_t t, size_t m, size_t p, Event e = Event::None, bool s=false, bool sNext=false, bool sNextNext=false) 
    : time(t), mode(m), phase(p), event(e), sto(s), stoNext(sNext), stoNextNext(sNextNext){};
};

/** Computes the time at which to interpolate, respecting interpolation rules around event times */
scalar_t getInterpolationTime(const Grid& grid);

/** Computes the interval start that respects interpolation rules around event times */
scalar_t getIntervalStart(const Grid& grid);

/** Computes the interval end that respects interpolation rules around event times */
scalar_t getIntervalEnd(const Grid& grid);

/** Computes the interval duration that respects interpolation rules around event times */
scalar_t getIntervalDuration(const Grid& start, const Grid& end);

std::vector<Grid> multiPhaseTimeDiscretizationGrid(scalar_t initTime, scalar_t finalTime, scalar_t dt, const ModeSchedule& modeSchedule,
                                                   const std::vector<bool>& isStoEnabled = {},
                                                   scalar_t dt_min = 10.0 * numeric_traits::limitEpsilon<scalar_t>());

/**
 * Decides on multi-phase time discretization along the horizon. Tries to makes step of dt, but will also ensure that eventtimes are part of the
 * discretization.
 *
 * @param initTime : start time.
 * @param finalTime : final time.
 * @param dt : desired discretization step.
 * @param eventTimes : Event times where a time discretization must be made.
 * @param dt_min : minimum discretization step. Smaller intervals will be merged. Needs to be bigger than limitEpsilon to avoid
 * interpolation problems
 * @return vector of discrete time points
 */
std::vector<AnnotatedTime> multiPhaseTimeDiscretization(scalar_t initTime, scalar_t finalTime, scalar_t dt,
                                                        const scalar_array_t& eventTimes,
                                                        scalar_t dt_min = 10.0 * numeric_traits::limitEpsilon<scalar_t>());

inline Grid::Event castEvent(const AnnotatedTime::Event& event) {
  switch (event)
  {
  case AnnotatedTime::Event::None:
    return Grid::Event::None;
    break;
  case AnnotatedTime::Event::PreEvent:
    return Grid::Event::PreEvent;
    break;
  case AnnotatedTime::Event::PostEvent:
    return Grid::Event::PostEvent;
    break;
  default:
    assert(false);
    break;
  }
}

std::ostream& operator<<(std::ostream& stream, const std::vector<AnnotatedTime>& timeDiscretization);

std::ostream& operator<<(std::ostream& stream, const std::vector<Grid>& timeDiscretizationGrid);

}  // namespace ocs2