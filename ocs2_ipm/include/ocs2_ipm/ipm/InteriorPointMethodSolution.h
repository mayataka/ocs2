#pragma once

#include <memory>

#include <ocs2_ipm/ipm/SlackDual.h>

namespace ocs2 {

/**
 * This class contains the solution related to the primal-dual interior point method.
 */
struct InteriorPointMethodSolution {
  /** Constructor */
  InteriorPointMethodSolution() = default;

  /** Destructor */
  ~InteriorPointMethodSolution() = default;

  /** Copy constructor */
  InteriorPointMethodSolution(const InteriorPointMethodSolution& other) = default;

  /** Copy Assignment */
  InteriorPointMethodSolution& operator=(const InteriorPointMethodSolution& other) = default;

  /** Move constructor */
  InteriorPointMethodSolution(InteriorPointMethodSolution&& other) noexcept = default;

  /** Move Assignment */
  InteriorPointMethodSolution& operator=(InteriorPointMethodSolution&& other) noexcept = default;

  /** Swap */
  void swap(InteriorPointMethodSolution& other) {
    slackDualTrajectory_.swap(other.slackDualTrajectory_);
  }

  void clear() {
    slackDualTrajectory_.clear();
  }

  std::vector<SlackDual> slackDualTrajectory_;
};


/**
 * Sets the default value to the given InteriorPointMethodSolution.
 *
 * @param[in] barrier: The barrier parameter. Must be positive.
 * @param[in] solution: Given InteriorPointMethodData.
 */
void setBarrier(scalar_t barrier, InteriorPointMethodSolution& solution);

/**
 * Checks that some coefficients of the given InteriorPointMethodSolution is positive.
 *
 * @param[in] solution: Given InteriorPointMethodSolution.
 * @param[in] solutionName: The name of the solution which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkPositive(const InteriorPointMethodSolution& solution, const std::string& solutionName);

}  // namespace ocs2
