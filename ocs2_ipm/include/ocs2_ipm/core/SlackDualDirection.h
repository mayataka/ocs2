#pragma once

#include <ocs2_core/Types.h>

namespace ocs2 {
namespace ipm {

/**
 * This class contains the directions of the slack and dual variables of the primal-dual interior point method.
 */
struct SlackDualDirection {
  /** Constructor */
  SlackDualDirection() = default;

  /** Destructor */
  ~SlackDualDirection() = default;

  /** Copy constructor */
  SlackDualDirection(const SlackDualDirection& other) = default;

  /** Copy Assignment */
  SlackDualDirection& operator=(const SlackDualDirection& other) = default;

  /** Move constructor */
  SlackDualDirection(SlackDualDirection&& other) noexcept = default;

  /** Move Assignment */
  SlackDualDirection& operator=(SlackDualDirection&& other) noexcept = default;

  /** Swap */
  void swap(SlackDualDirection& other) {
    slackDirection.swap(other.slackDirection);
    dualDirection.swap(other.dualDirection);
  }

  void resize(size_t size) {
    slackDirection.resize(size);
    dualDirection.resize(size);
  }

  vector_t slackDirection;
  vector_t dualDirection;
};

}  // namespace ipm
}  // namespace ocs2
