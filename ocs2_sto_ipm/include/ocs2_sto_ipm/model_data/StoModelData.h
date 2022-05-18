#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "ocs2_core/Types.h"
#include "ocs2_sto_ipm/Types.h"

namespace ocs2 {
namespace sto_ipm {

/**
 * The optimal control problem with switching time optimization model data.
 */
struct StoModelData {
  int stateDim = 0;
  int inputDim = 0;
  scalar_t time = 0.0;

  // Dynamics
  vector_t dynamicsBias;
  matrix_t dynamicsCovariance;
  VectorFunctionLinearApproximationWrapper dynamics;

  // Cost
  ScalarFunctionQuadraticApproximationWrapper cost;

  // Inequality constraints
  VectorFunctionLinearApproximationWrapper stateIneqConstraint;
  VectorFunctionLinearApproximationWrapper inputIneqConstraint;
  VectorFunctionLinearApproximationWrapper stateInputIneqConstraint;

  // Equality constraints
  VectorFunctionLinearApproximationWrapper stateEqConstraint;
  VectorFunctionLinearApproximationWrapper stateInputEqConstraint;
};

/**
 * Checks whether the size of the StoModelData variables matches the given state and input dimensions.
 *
 * @param [in] data: The StoModelData to be examined.
 * @param [in] stateDim: The dimension of the state vector.
 * @param [in] inputDim: The dimension of the input vector.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string checkSize(const StoModelData& data, int stateDim, int inputDim);

/**
 * Checks the numerical properties of the cost function and its derivatives.
 *
 * @param [in] data: The StoModelData to be examined.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string checkCostProperties(const StoModelData& data);

/**
 * Checks if the Shur complement of cost Hessian w.r.t. input is psd.
 *
 * @param [in] cost: Cost function quadratic approximation.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string schurComplementOfCostHessianIsPsd(const ScalarFunctionQuadraticApproximationWrapper& cost);

/**
 * Checks the numerical properties of the dynamics derivatives.
 *
 * @param [in] data: The StoModelData to be examined.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string checkDynamicsProperties(const StoModelData& data);

/**
 * Checks if the linearized system is controllable.
 *
 * @param [in] dynamics: Dynamics linear approximation.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string checkControllability(const VectorFunctionLinearApproximationWrapper& dynamics);

/**
 * Checks the numerical properties of the constraint functions and derivatives.
 *
 * @param [in] data: The StoModelData to be examined.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string checkConstraintProperties(const StoModelData& data);

}  // namespace sto_ipm
}  // namespace ocs2
