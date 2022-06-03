#pragma once

#include <memory>

#include <ocs2_core/Types.h>

namespace ocs2 {

/**
 * This class contains the slack and dual variables of the primal-dual interior point method.
 */
struct SlackDual {
  /** Constructor */
  SlackDual() = default;

  /** Destructor */
  ~SlackDual() = default;

  /** Copy constructor */
  SlackDual(const SlackDual& other) = default;

  /** Copy Assignment */
  SlackDual& operator=(const SlackDual& other) = default;

  /** Move constructor */
  SlackDual(SlackDual&& other) noexcept = default;

  /** Move Assignment */
  SlackDual& operator=(SlackDual&& other) noexcept = default;

  /** Swap */
  void swap(SlackDual& other) {
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
 * @param[in] slackDual: Given SlackDual.
 */
void setBarrier(size_t nc, scalar_t barrier, SlackDual& slackDual);

/**
 * Sets the default values to slack and dual.
 *
 * @param[in] barrier: The barrier parameter. Must be positive.
 * @param[in] slackDual: Given SlackDual.
 */
void setBarrier(scalar_t barrier, SlackDual& slackDual);

/**
 * Checks the size of the given SlackDual.
 *
 * @param[in] constraintDim: Constraint dimension.
 * @param[in] slackDual: Given SlackDual.
 * @param[in] dataName: The name of the data which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkSize(int constraintDim, const SlackDual& slackDual, const std::string& dataName);

/**
 * Checks that some coefficients of the given SlackDual is positive.
 *
 * @param[in] slackDual: Given SlackDual.
 * @param[in] dataName: The name of the data which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkPositive(const SlackDual& slackDual, const std::string& dataName);

}  // namespace ocs2
