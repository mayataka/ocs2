#include <ocs2_ipm/approximate_model/ConstraintProjection.h>
#include <ocs2_sqp/ConstraintProjection.h>
#include <ocs2_oc/approximate_model/ChangeOfInputVariables.h>
#include <ocs2_ipm/approximate_model/ChangeOfInputVariables.h>

namespace ocs2 {
namespace ipm {

void projectIntermediateLQ(const ModelData& modelData, VectorFunctionLinearApproximation& constraintProjection, 
                           ModelData& projectedModelData) {
  const auto constraintProjectionDim = modelData.stateInputEqConstraint.f.size();

  // dimensions and time
  projectedModelData.stateDim = modelData.stateDim;
  projectedModelData.inputDim = modelData.inputDim - constraintProjectionDim;
  projectedModelData.time     = modelData.time;

  // dynamics
  projectedModelData.dynamicsBias = modelData.dynamicsBias;
  projectedModelData.dynamicsCovariance = modelData.dynamicsCovariance;
  projectedModelData.dynamics = modelData.dynamics;

  // cost
  projectedModelData.cost = modelData.cost;

  // Equality constraints
  projectedModelData.stateInputEqConstraint = modelData.stateInputEqConstraint;
  projectedModelData.stateEqConstraint = modelData.stateEqConstraint;

  // inequality constraints
  projectedModelData.stateIneqConstraint = modelData.stateIneqConstraint;
  projectedModelData.stateInputIneqConstraint = modelData.stateInputIneqConstraint;

  // hamiltonian
  projectedModelData.hamiltonian = modelData.hamiltonian;

  if (constraintProjectionDim > 0) {
    // Projection stored instead of constraint, 
    constraintProjection = luConstraintProjection(modelData.stateInputEqConstraint);
    // Adapt dynamics and cost
    changeOfInputVariables(projectedModelData.dynamics, constraintProjection.dfdu, constraintProjection.dfdx, constraintProjection.f);
    changeOfInputVariables(projectedModelData.cost, constraintProjection.dfdu, constraintProjection.dfdx, constraintProjection.f);
    changeOfInputVariables(projectedModelData.hamiltonian, constraintProjection.dfdu, constraintProjection.dfdx, constraintProjection.f);
  }
}

}  // namespace ipm
}  // namespace ocs2