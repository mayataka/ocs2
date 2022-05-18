#pragma once

#include <ostream>
#include <vector>

#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>
#include <ocs2_sto_ipm/TimeDiscretization.h>


namespace ocs2 {
namespace sto_ipm {

/**
 * Defines modes at each discretized grid.
 */
struct DiscreteTimeModeSchedule {
  /**
   * Default constructor.
   */
  DiscreteTimeModeSchedule() : DiscreteTimeModeSchedule(std::vector<size_t>{0}, std::vector<bool>{false}) {}

  /**
   * Constructor for a DiscreteTimeModeSchedule. The number of modes must be greater than zero (N > 0)
   * @param [in] modeSequenceInput : mode sequence of size N.
   * @param [in] isStoEnabled : isStoEnabled of size modeSequenceInput.back()+1.
   */
  DiscreteTimeModeSchedule(std::vector<size_t> modeSequenceInput, std::vector<bool> isStoEnabled);

  /**
   *  Returns the mode of the specified time stage.
   *  @param [in] timeStage: The inquiry time stage.
   *  @return the associated mode for the input time stage.
   */
  size_t modeAtTimeStage(size_t timeStage) const;

  /**
   *  Returns the phase of the specified time stage.
   *  @param [in] timeStage: The inquiry time stage.
   *  @return the associated phase for the input time stage.
   */
  size_t phaseAtTimeStage(size_t timeStage) const;

  /**
   *  Returns if the switching time optimization (STO) is enabled at the specified phase.
   *  The phase is counted from the beginning of the mode sequence.
   *  If the input phase is larger than the total number of phases, return false.
   *  @param [in] phase: The inquiry phase.
   *  @return if the switching time optimization (STO) is enabled or not.
   */
  bool isStoEnabledAtPhase(size_t phase) const;

  /** Clears modeSchedule */
  void clear() {
    modeSequence.clear();
    phaseSequence.clear();
    isStoEnabled.clear();
  }

  std::vector<size_t> modeSequence;  // mode sequence of size N+1 (N: length of the horizon)
  std::vector<size_t> phaseSequence;  // mode index sequence of size N+1 (N: length of the horizon)
  std::vector<bool> isStoEnabled;     // sequence of STO flags of size N+1
};

/** Exchanges the given values. */
void swap(DiscreteTimeModeSchedule& lh, DiscreteTimeModeSchedule& rh);

/** Inserts modeSchedule into the output stream. */
std::ostream& operator<<(std::ostream& stream, const DiscreteTimeModeSchedule& modeSchedule);

}  // namespace sto_ipm
}  // namespace ocs2
