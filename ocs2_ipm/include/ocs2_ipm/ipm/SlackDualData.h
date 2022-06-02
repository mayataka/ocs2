#pragma once

#include <memory>

#include <ocs2_core/Types.h>

namespace ocs2 {

/**
 * This class contains the slack and dual variables of the primal-dual interior point method.
 */
struct SlackDualData {
  /** Constructor */
  SlackDualData() = default;

  /** Destructor */
  ~SlackDualData() = default;

  /** Copy constructor */
  SlackDualData(const SlackDualData& other) = default;

  /** Copy Assignment */
  SlackDualData& operator=(const SlackDualData& other) = default;

  /** Move constructor */
  SlackDualData(SlackDualData&& other) noexcept = default;

  /** Move Assignment */
  SlackDualData& operator=(SlackDualData&& other) noexcept = default;

  /** Swap */
  void swap(SlackDualData& other) {
    slack.swap(other.slack);
    dual.swap(other.dual);
  }

  scalar_t barrier;
  vector_t slack;
  vector_t dual;
};


/**
 * Resizes and sets the default values to slack and dual.
 *
 * @param[in] nc: Size of the constraint. Must be non-negative.
 * @param[in] barrier: The barrier parameter. Must be positive.
 * @param[in] slackDual: Given SlackDualData.
 */
void setBarrier(size_t nc, scalar_t barrier, SlackDualData& slackDual);

/**
 * Sets the default values to slack and dual.
 *
 * @param[in] barrier: The barrier parameter. Must be positive.
 * @param[in] slackDual: Given SlackDualData.
 */
void setBarrier(scalar_t barrier, SlackDualData& slackDual);

/**
 * Checks the size of the given SlackDualData.
 *
 * @param[in] constraintDim: Constraint dimension.
 * @param[in] slackDual: Given SlackDualData.
 * @param[in] dataName: The name of the data which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkSize(int constraintDim, const SlackDualData& slackDual, const std::string& dataName);

/**
 * Checks that some coefficients of the given SlackDualData is positive.
 *
 * @param[in] slackDual: Given SlackDualData.
 * @param[in] dataName: The name of the data which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkPositive(const SlackDualData& slackDual, const std::string& dataName);

}  // namespace ocs2
