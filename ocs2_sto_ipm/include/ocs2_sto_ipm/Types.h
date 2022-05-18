#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <ostream>
#include <vector>
#include <ocs2_core/Types.h>

namespace ocs2 {
namespace sto_ipm {

/**
 * Defines the quadratic approximation
 * f(x,u) = 1/2 dx' dfdxx dx + du' dfdux dx + 1/2 du' dfduu du + dfdx' dx + dfdu' du + f
 * Also includes the sensitivities of the discretized function w.r.t the switching times.
 */
struct ScalarFunctionQuadraticApproximationWrapper {
  /** Base quadratic approximation w.r.t. state and input */
  ScalarFunctionQuadraticApproximation base;
  /** Second derivative w.r.t state (lhs) and time (rhs) */
  vector_t dfdxt;
  /** Second derivative w.r.t input (lhs) and time (rhs) */
  vector_t dfdut;
  /** Second derivative w.r.t time */
  scalar_t dfdtt;
  /** First derivative w.r.t time */
  scalar_t dfdt;

  /** Default constructor */
  ScalarFunctionQuadraticApproximationWrapper() = default;

  /** Construct and resize the members to given size. */
  ScalarFunctionQuadraticApproximationWrapper(size_t nx, size_t nu);

  /** Compound addition assignment operator */
  ScalarFunctionQuadraticApproximationWrapper& operator+=(const ScalarFunctionQuadraticApproximationWrapper& rhs);

  /** Compound scalar multiplication and assignment operator */
  ScalarFunctionQuadraticApproximationWrapper& operator*=(scalar_t scalar);

  /**
   * Resize the members to the given size
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   */
  ScalarFunctionQuadraticApproximationWrapper& resize(size_t nx, size_t nu);

  /**
   * Resizes the members to the given size, and sets all coefficients to zero.
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   */
  ScalarFunctionQuadraticApproximationWrapper& setZero(size_t nx, size_t nu);

  /**
   * Factory function with zero initialization
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   * @return Zero initialized object of given size.
   */
  static ScalarFunctionQuadraticApproximationWrapper Zero(size_t nx, size_t nu);
};

std::ostream& operator<<(std::ostream& out, const ScalarFunctionQuadraticApproximationWrapper& f);

/**
 * Checks that the given quadratic approximation is valid, self-adjoint, and positive semi-definite (PSD).
 * @param[in] data: Given quadratic approximation.
 * @param[in] dataName: The name of the data which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkBeingPSD(const ScalarFunctionQuadraticApproximationWrapper& data, const std::string& dataName);

/**
 * Checks the size of the given quadratic approximation.
 *
 * @param[in] stateDim: Number of states.
 * @param[in] inputDim: Number of inputs.
 * @param[in] data: Given quadratic approximation.
 * @param[in] dataName: The name of the data which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkSize(int stateDim, int inputDim, const ScalarFunctionQuadraticApproximationWrapper& data, const std::string& dataName);

inline ScalarFunctionQuadraticApproximationWrapper operator*(ScalarFunctionQuadraticApproximationWrapper lhs, scalar_t scalar) {
  return lhs *= scalar;
}
inline ScalarFunctionQuadraticApproximationWrapper operator*(scalar_t scalar, ScalarFunctionQuadraticApproximationWrapper rhs) {
  return rhs *= scalar;
}

/**
 * Defines the linear model of a vector-valued function
 * f(x,u) = dfdx dx + dfdu du + f
 * Also includes the sensitivities of the discretized function w.r.t the switching times.
 */
struct VectorFunctionLinearApproximationWrapper {
  /** Base linear approximation w.r.t. state and input */
  VectorFunctionLinearApproximation base;
  /** Derivative w.r.t time */
  vector_t dfdt;

  /** Default constructor */
  VectorFunctionLinearApproximationWrapper() = default;

  /** Construct and resize the members to given size. */
  VectorFunctionLinearApproximationWrapper(size_t nv, size_t nx, size_t nu);

  /**
   * Resize the members to the given size
   * @param[in] nv Vector dimension
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   */
  VectorFunctionLinearApproximationWrapper& resize(size_t nv, size_t nx, size_t nu);

  /**
   * Resizes the members to the given size, and sets all coefficients to zero.
   * @param[in] nv Vector dimension
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   */
  VectorFunctionLinearApproximationWrapper& setZero(size_t nv, size_t nx, size_t nu);

  /**
   * Factory function with zero initialization
   * @param[in] nv Vector dimension
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   * @return Zero initialized object of given size.
   */
  static VectorFunctionLinearApproximationWrapper Zero(size_t nv, size_t nx, size_t nu);
};

std::ostream& operator<<(std::ostream& out, const VectorFunctionLinearApproximationWrapper& f);

/**
 * Checks the size of the given vector-function linear approximation.
 *
 * @param[in] vectorDim: The vector function dimension.
 * @param[in] stateDim: Number of states.
 * @param[in] inputDim: Number of inputs.
 * @param[in] data: Given linear approximation.
 * @param[in] dataName: The name of the data which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkSize(int vectorDim, int stateDim, int inputDim, const VectorFunctionLinearApproximationWrapper& data,
                      const std::string& dataName);

/**
 * Defines quadratic approximation of a vector-valued function
 * f[i](x,u) = 1/2 dx' dfdxx[i] dx + du' dfdux[i] dx + 1/2 du' dfduu[i] du + dfdx[i,:] dx + dfdu[i,:] du + f[i]
 * Also includes the sensitivities of the discretized function w.r.t the switching times.
 */
struct VectorFunctionQuadraticApproximationWrapper {
  /** Base quadratic approximation w.r.t. state and input */
  VectorFunctionQuadraticApproximation base;
  /** Second derivative w.r.t state (lhs) and time (rhs) */
  matrix_t dfdxt;
  /** Second derivative w.r.t input (lhs) and time (rhs) */
  matrix_t dfdut;
  /** Second derivative w.r.t time */
  vector_t dfdtt;
  /** First derivative w.r.t time */
  vector_t dfdt;

  /** Default constructor */
  VectorFunctionQuadraticApproximationWrapper() = default;

  /** Construct and resize the members to given size. */
  VectorFunctionQuadraticApproximationWrapper(size_t nv, size_t nx, size_t nu);

  /**
   * Resize the members to the given size
   * @param[in] nv Vector dimension
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   */
  VectorFunctionQuadraticApproximationWrapper& resize(size_t nv, size_t nx, size_t nu);

  /**
   * Resizes the members to the given size, and sets all coefficients to zero.
   * @param[in] nv Vector dimension
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   */
  VectorFunctionQuadraticApproximationWrapper& setZero(size_t nv, size_t nx, size_t nu);

  /**
   * Factory function with zero initialization
   * @param[in] nv Vector dimension
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   * @return Zero initialized object of given size.
   */
  static VectorFunctionQuadraticApproximationWrapper Zero(size_t nv, size_t nx, size_t nu);
};

std::ostream& operator<<(std::ostream& out, const VectorFunctionQuadraticApproximationWrapper& f);

}  // namespace sto_ipm
}  // namespace ocs2
