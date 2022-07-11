#pragma once

#include <string>

#include <ocs2_core/Types.h>

namespace ocs2 {
namespace ipm {

/** Data for the interior point method */
struct InteriorPointMethodData {
  // Dimension of the constraint
  int dim;

  // Primal residual and complementary slackness
  vector_t primalResidual;
  vector_t complementary;

  // Slack and dual search directions
  vector_t slackDirection;
  vector_t dualDirection;

  // Temporal data in condensing
  vector_t cond, dualDivSlack;
  VectorFunctionLinearApproximation linearApproximation;

  /** Default constructor */
  InteriorPointMethodData() = default;

  /** Construct and resize the members to given size. */
  InteriorPointMethodData(size_t nc, size_t nx, size_t nu);

  /**
   * Resize the members to the given size
   * @param[in] nc Constraint dimension
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   */
  InteriorPointMethodData& resize(size_t nc, size_t nx, size_t nu);

  /**
   * Resizes the members to the given size, and sets all coefficients to zero
   * @param[in] nc Constraint dimension
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   */
  InteriorPointMethodData& setZero(size_t nc, size_t nx, size_t nu);

  /**
   * Factory function with zero initialization
   * @param[in] nc Constraint dimension
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   * @return Zero initialized object of given size.
   */
  static InteriorPointMethodData Zero(size_t nc, size_t nx, size_t nu);
};

/**
 * Checks the size of the given InteriorPointMethodData.
 *
 * @param[in] constraintDim: Constraint dimension.
 * @param[in] data: Given InteriorPointMethodData.
 * @param[in] dataName: The name of the data which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkSize(int constraintDim, const InteriorPointMethodData& data, const std::string& dataName);

}  // namespace ipm
}  // namespace ocs2
