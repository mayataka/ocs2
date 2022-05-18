#pragma once

#include <memory>

#include <ocs2_core/Types.h>
#include <ocs2_sto_ipm/Types.h>

namespace ocs2 {
namespace sto_ipm {

/** Data for the interior point method */
struct InteriorPointMethodData {
  // Dimension of the constraint
  int dim;

  // The barrier parameter
  scalar_t barrier;

  // Slack and dual variables
  vector_t slack;
  vector_t dual;

  // Primal residual and complementary slackness
  vector_t primalResidual;
  vector_t complementary;

  // Slack and dual search directions
  vector_t slackDirection;
  vector_t dualDirection;

  // Temporal data in condensing
  vector_t cond;
  VectorFunctionLinearApproximation linearApproximation;

  /** Default constructor */
  InteriorPointMethodData() = default;

  /** Construct and resize the members to given size. */
  InteriorPointMethodData(size_t nc, size_t nx, size_t nu);

  /** Construct and resize the members to given size. Also set the barrier parameter */
  InteriorPointMethodData(size_t nc, size_t nx, size_t nu, scalar_t barrier);

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
   * Resizes the members to the given size, and sets some coefficients to using
   * the barrier parameter and the others by zero.
   * @param[in] barrier Barrier parameter. 
   */
  InteriorPointMethodData& setBarrier(scalar_t barrier);

  /**
   * Resizes the members to the given size, and sets some coefficients to using
   * the barrier parameter and the others by zero.
   * @param[in] nc Constraint dimension
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   * @param[in] barrier Barrier parameter. 
   */
  InteriorPointMethodData& setBarrier(size_t nc, size_t nx, size_t nu, scalar_t barrier);

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

/**
 * Checks that some coefficients of the given InteriorPointMethodData is positive.
 *
 * @param[in] data: Given InteriorPointMethodData.
 * @param[in] dataName: The name of the data which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkPositive(const InteriorPointMethodData& data, const std::string& dataName);

}  // namespace sto_ipm
}  // namespace ocs2
