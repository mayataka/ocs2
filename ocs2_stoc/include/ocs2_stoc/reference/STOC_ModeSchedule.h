#pragma once

#include <ostream>
#include <vector>

#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>

namespace ocs2 {
namespace stoc {

/**
 * Defines a sequence of N modes, separated by N-1 event times
 */
struct STOC_ModeSchedule {
  /**
   * Default constructor.
   */
  STOC_ModeSchedule() : STOC_ModeSchedule(std::vector<scalar_t>{}, std::vector<size_t>{0}, std::vector<bool>{}) {}

  /**
   * Constructor for a ModeSchedule. The number of modes must be greater than zero (N > 0)
   * @param [in] eventTimesInput : event times of size N - 1
   * @param [in] modeSequenceInput : mode sequence of size N
   * @param [in] isStoEnabledInput: STO sequence of size N - 1
   */
  STOC_ModeSchedule(std::vector<scalar_t> eventTimesInput, std::vector<size_t> modeSequenceInput, std::vector<bool> isStoEnabledInput);

  /**
   *  Returns the mode based on the query time.
   *  Events are counted as follows:
   *      ------ | ------ | ------ | ...  ------ | ------
   *         t[0]     t[1]     t[2]        t[n-1]
   *  mode: m[0]    m[1]     m[2] ...     m[n-1]    m[n]
   *
   *  If time equal to a switch time is requested, the lower count is taken
   *
   *  @param [in] time: The inquiry time.
   *  @return the associated mode for the input time.
   */
  size_t modeAtTime(scalar_t time) const;

  /** Clears modeSchedule */
  void clear() {
    eventTimes.clear();
    modeSequence.clear();
    isStoEnabled.clear();
  }

  std::vector<scalar_t> eventTimes;  // event times of size N - 1
  std::vector<size_t> modeSequence;  // mode sequence of size N
  std::vector<bool> isStoEnabled;  // STO sequence of size N - 1
};

/** Exchanges the given values. */
void swap(STOC_ModeSchedule& lh, STOC_ModeSchedule& rh);

/** Inserts modeSchedule into the output stream. */
std::ostream& operator<<(std::ostream& stream, const STOC_ModeSchedule& modeSchedule);

}  // namespace stoc
}  // namespace ocs2
