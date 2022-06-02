#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <ocs2_core/Types.h>
#include <ocs2_core/model_data/ModelData.h>
#include <ocs2_ipm/model_data/Hamiltonian.h>
#include <ocs2_ipm/ipm/InteriorPointMethodData.h>

namespace ocs2 {
namespace ipm {

/**
 * The optimal control problem model data for interior point method.
 */
struct ModelData {
  int stateDim = 0;
  int inputDim = 0;
  scalar_t time = 0.0;

  // Dynamics
  vector_t dynamicsBias;
  matrix_t dynamicsCovariance;
  VectorFunctionLinearApproximation dynamics;

  // Cost
  ScalarFunctionQuadraticApproximation cost;

  // Equality constraints
  VectorFunctionLinearApproximation stateEqConstraint;
  VectorFunctionLinearApproximation stateInputEqConstraint;

  // Inequality constraints
  VectorFunctionLinearApproximation stateIneqConstraint;
  VectorFunctionLinearApproximation stateInputIneqConstraint;

  // Interior point method data
  InteriorPointMethodData stateIneqIpmData;
  InteriorPointMethodData stateInputIneqIpmData;

  // Hamiltonian
  Hamiltonian hamiltonian;
};

/**
 * Checks whether the size of the ModelData variables matches the given state and input dimensions.
 *
 * @param [in] data: The ModelData to be examined.
 * @param [in] stateDim: The dimension of the state vector.
 * @param [in] inputDim: The dimension of the input vector.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string checkSize(const ModelData& data, int stateDim, int inputDim);

/**
 * Checks the numerical properties of the cost function and its derivatives.
 *
 * @param [in] data: The ModelData to be examined.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string checkCostProperties(const ModelData& data);

/**
 * Checks if the Shur complement of cost Hessian w.r.t. input is psd.
 *
 * @param [in] cost: Cost function quadratic approximation.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string schurComplementOfCostHessianIsPsd(const ScalarFunctionQuadraticApproximation& cost);

/**
 * Checks the numerical properties of the dynamics derivatives.
 *
 * @param [in] data: The ModelData to be examined.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string checkDynamicsProperties(const ModelData& data);

/**
 * Checks if the linearized system is controllable.
 *
 * @param [in] dynamics: Dynamics linear approximation.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string checkControllability(const VectorFunctionLinearApproximation& dynamics);

/**
 * Checks the numerical properties of the constraint functions and derivatives.
 *
 * @param [in] data: The ModelData to be examined.
 * @return The description of the error. If there was no error it would be empty.
 */
std::string checkConstraintProperties(const ModelData& data);

}  // namespace ipm
}  // namespace ocs2
