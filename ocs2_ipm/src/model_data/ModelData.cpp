#include "ocs2_ipm/model_data/ModelData.h"

#include <iostream>

namespace ocs2 {
namespace ipm {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkSize(const ModelData& data, int stateDim, int inputDim) {
  std::stringstream errorDescription;

  errorDescription << checkSize(static_cast<ModelData>(data), stateDim, inputDim);

  // hamiltonian
  errorDescription << checkSize(stateDim, inputDim, data.hamiltonian, "hamiltonian");

  // state inequality constraints
  errorDescription << checkSize(data.stateIneqConstraint.f.size(), stateDim, 0, data.stateIneqConstraint, "stateIneqConstraint");

  // state-input inequality constraints
  errorDescription << checkSize(data.stateInputIneqConstraint.f.size(), stateDim, inputDim, data.stateInputIneqConstraint, 
                                "stateInputIneqConstraint");

  return errorDescription.str();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkCostProperties(const ModelData& data) {
  std::stringstream errorDescription;

  errorDescription << checkCostProperties(static_cast<ModelData>(data));

  return errorDescription.str();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkDynamicsProperties(const ModelData& data) {
  std::stringstream errorDescription;

  errorDescription << checkDynamicsProperties(static_cast<ModelData>(data));

  return errorDescription.str();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkConstraintProperties(const ModelData& data) {
  std::stringstream errorDescription;

  errorDescription << checkConstraintProperties(static_cast<ModelData>(data));

  if (data.stateIneqConstraint.f.rows() > 0) {
    if (!data.stateIneqConstraint.f.allFinite()) {
      errorDescription << "State-only inequality constraint is not finite.\n";
    }
    if (!data.stateIneqConstraint.dfdx.allFinite()) {
      errorDescription << "State-only inequality constraint derivative w.r.t. state is not finite.\n";
    }
  }

  if (data.stateInputIneqConstraint.f.rows() > 0) {
    const auto inputDim = data.stateInputIneqConstraint.dfdu.cols();
    const auto numConstraints = data.stateInputIneqConstraint.f.rows();

    if (!data.stateInputIneqConstraint.f.allFinite()) {
      errorDescription << "Input-state inequality constraint is not finite.\n";
    }
    if (!data.stateInputIneqConstraint.dfdx.allFinite()) {
      errorDescription << "Input-state inequality constraint derivative w.r.t. state is not finite.\n";
    }
    if (!data.stateInputIneqConstraint.dfdu.allFinite()) {
      errorDescription << "Input-state inequality constraint derivative w.r.t. input is not finite.\n";
    }
  }

  return errorDescription.str();
}


std::string checkHamiltonianProperties(const ModelData& data) {
  std::stringstream errorDescription;
  
  return errorDescription.str();
}

}  // namespace ipm
}  // namespace ocs2
