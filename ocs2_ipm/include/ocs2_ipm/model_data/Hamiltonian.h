#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <ocs2_core/Types.h>

namespace ocs2 {

struct Hamiltonian {
  /** Value of the hamiltonian. */
  scalar_t h = 0.;
  /** Derivative of the hamiltonian w.r.t. the switching time. */
  scalar_t dhdt;
  /** Derivative of the hamiltonian w.r.t. the state. */
  vector_t dhdx;
  /** Derivative of the hamiltonian w.r.t. the input. */
  vector_t dhdu;
  /** The derivative of the state equation w.r.t. the switching time. 
   * That is, the derivative of the hamiltonian w.r.t. the Lagrange multiplier. */
  vector_t dfdt;

  /** Default constructor */
  Hamiltonian() = default;

  /** Construct and resize the members to given size. */
  Hamiltonian(size_t nx, size_t nu);

  /** Compound addition assignment operator */
  Hamiltonian& operator+=(const Hamiltonian& rhs);

  /** Compound scalar multiplication and assignment operator */
  Hamiltonian& operator*=(scalar_t scalar);

  /**
   * Resize the members to the given size
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   */
  Hamiltonian& resize(size_t nx, size_t nu);

  /**
   * Resizes the members to the given size, and sets all coefficients to zero.
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   */
  Hamiltonian& setZero(size_t nx, size_t nu);

  /**
   * Factory function with zero initialization
   * @param[in] nx State dimension
   * @param[in] nu Input dimension
   * @return Zero initialized object of given size.
   */
  static Hamiltonian Zero(size_t nx, size_t nu);
};

std::ostream& operator<<(std::ostream& out, const Hamiltonian& hamiltonian);

/**
 * Checks the size of the given quadratic approximation.
 *
 * @param[in] stateDim: Number of states.
 * @param[in] inputDim: Number of inputs.
 * @param[in] hamiltonian: Given hamiltonian data.
 * @param[in] hamiltonianName: The name of the hamiltonian data which appears in the output error message.
 * @return The description of the error. If there was no error it would be empty;
 */
std::string checkSize(int stateDim, int inputDim, const Hamiltonian& data, const std::string& dataName);

inline Hamiltonian operator*(Hamiltonian lhs, scalar_t scalar) {
  return lhs *= scalar;
}

inline Hamiltonian operator*(scalar_t scalar, Hamiltonian rhs) {
  return rhs *= scalar;
}

}  // namespace ocs2
