#include <iostream>

#include "ocs2_core/misc/LinearAlgebra.h"
#include "ocs2_sto_ipm/model_data/StoModelData.h"

namespace ocs2 {
namespace sto_ipm {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkSize(const StoModelData& data, int stateDim, int inputDim) {
  std::stringstream errorDescription;

  if (data.stateDim != stateDim) {
    errorDescription << "data.stateDim != " << stateDim << "\n";
  }
  if (data.inputDim != inputDim) {
    errorDescription << "data.inputDim != " << inputDim << "\n";
  }

  // dynamics
  if (data.dynamics.base.f.size() > 0) {
    errorDescription << checkSize(stateDim, stateDim, inputDim, data.dynamics, "dynamics");

    if (data.dynamicsBias.size() != stateDim) {
      errorDescription << "dynamicsBias.size() != " << stateDim << "\n";
    }
  }

  // cost
  errorDescription << checkSize(stateDim, inputDim, data.cost, "cost");

  // state inequality constraints
  errorDescription << checkSize(data.stateIneqConstraint.base.f.size(), stateDim, 0, data.stateIneqConstraint, "stateIneqConstraint");

  // state inequality constraints
  errorDescription << checkSize(data.stateInputIneqConstraint.base.f.size(), stateDim, inputDim, data.stateInputIneqConstraint, 
                                "stateInputIneqConstraint");

  // state equality constraints
  errorDescription << checkSize(data.stateEqConstraint.base.f.size(), stateDim, 0, data.stateEqConstraint, "stateEqConstraint");

  // state-input equality constraints
  errorDescription << checkSize(data.stateInputEqConstraint.base.f.size(), stateDim, inputDim, data.stateInputEqConstraint,
                                "stateInputEqConstraint");

  return errorDescription.str();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkCostProperties(const StoModelData& data) {
  std::stringstream errorDescription;

  errorDescription << checkBeingPSD(data.cost, "cost");

  // check if R is invertible and its schur complement is PSD
  if (data.cost.base.dfduu.size() > 0) {
    const auto rcond = data.cost.base.dfduu.ldlt().rcond();
    if (rcond < Eigen::NumTraits<scalar_t>::epsilon()) {
      errorDescription << "Cost second derivative w.r.t. input is not invertible. It's reciprocal condition number is " +
                              std::to_string(rcond) + ".\n";
    } else {
      // check schur complement of R being PSD
      errorDescription << schurComplementOfCostHessianIsPsd(data.cost);
    }
  }

  return errorDescription.str();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string schurComplementOfCostHessianIsPsd(const ScalarFunctionQuadraticApproximationWrapper& cost) {
  if (cost.base.dfdxx.size() > 0 && cost.base.dfduu.size() > 0) {
    matrix_t UofUUT;
    LinearAlgebra::computeInverseMatrixUUT(cost.base.dfduu, UofUUT);
    const matrix_t UT_P = UofUUT.transpose() * cost.base.dfdux;
    matrix_t inputHessianSchurComplement = cost.base.dfdxx;
    inputHessianSchurComplement.noalias() -= UT_P.transpose() * UT_P;

    // check for being psd
    return ::ocs2::checkBeingPSD(inputHessianSchurComplement, "Schur complement of cost second derivative w.r.t. input");

  } else {
    return "Either cost.dfdxx or cost.dfduu are not set!";
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkDynamicsProperties(const StoModelData& data) {
  std::stringstream errorDescription;

  if (!data.dynamics.base.f.allFinite()) {
    errorDescription << "Dynamics is not finite.";
  }
  if (!data.dynamicsBias.allFinite()) {
    errorDescription << "Dynamics bias is not finite.";
  }
  if (!data.dynamics.base.dfdx.allFinite()) {
    errorDescription << "Dynamics derivative w.r.t. state is not finite.";
  }
  if (!data.dynamics.base.dfdu.allFinite()) {
    errorDescription << "Dynamics derivative w.r.t. input is not finite.";
  }
  if (!data.dynamics.dfdt.allFinite()) {
    errorDescription << "Dynamics derivative w.r.t. switching time is not finite.";
  }

  return errorDescription.str();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkControllability(const VectorFunctionLinearApproximationWrapper& dynamics) {
  const size_t stateDim = dynamics.base.dfdu.rows();
  const size_t inputDim = dynamics.base.dfdu.cols();

  // controllability matrix
  matrix_t ctrlMatrix(stateDim, inputDim * stateDim);
  ctrlMatrix.leftCols(inputDim) = dynamics.base.dfdu;
  for (size_t i = 1; i < stateDim; i++) {
    ctrlMatrix.middleCols(i * inputDim, inputDim).noalias() = dynamics.base.dfdx * ctrlMatrix.middleCols((i - 1) * inputDim, inputDim);
  }

  // controllability rank
  const size_t ctrlMatrixRank = LinearAlgebra::rank(ctrlMatrix);

  std::stringstream errorDescription;
  if (ctrlMatrixRank < stateDim) {
    errorDescription << "Uncontrollable system: controllability matrix rank should be " << stateDim << " as opposed to " << ctrlMatrixRank
                     << "\n";
  }

  return errorDescription.str();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkConstraintProperties(const StoModelData& data) {
  std::stringstream errorDescription;

  if (data.stateIneqConstraint.base.f.rows() > 0) {
    if (!data.stateIneqConstraint.base.f.allFinite()) {
      errorDescription << "State-only inequality constraint is not finite.\n";
    }
    if (!data.stateIneqConstraint.base.dfdx.allFinite()) {
      errorDescription << "State-only inequality constraint derivative w.r.t. state is not finite.\n";
    }
  }

  if (data.stateInputIneqConstraint.base.f.rows() > 0) {
    const auto inputDim = data.stateInputIneqConstraint.base.dfdu.cols();
    const auto numConstraints = data.stateInputIneqConstraint.base.f.rows();

    if (!data.stateInputIneqConstraint.base.f.allFinite()) {
      errorDescription << "Input-state inequality constraint is not finite.\n";
    }
    if (!data.stateInputIneqConstraint.base.dfdx.allFinite()) {
      errorDescription << "Input-state inequality constraint derivative w.r.t. state is not finite.\n";
    }
    if (!data.stateInputIneqConstraint.base.dfdu.allFinite()) {
      errorDescription << "Input-state inequality constraint derivative w.r.t. input is not finite.\n";
    }
  }

  if (data.stateEqConstraint.base.f.rows() > 0) {
    if (!data.stateEqConstraint.base.f.allFinite()) {
      errorDescription << "State-only equality constraint is not finite.\n";
    }
    if (!data.stateEqConstraint.base.dfdx.allFinite()) {
      errorDescription << "State-only equality constraint derivative w.r.t. state is not finite.\n";
    }
  }

  if (data.stateInputEqConstraint.base.f.rows() > 0) {
    const auto inputDim = data.stateInputEqConstraint.base.dfdu.cols();
    const auto numConstraints = data.stateInputEqConstraint.base.f.rows();

    if (!data.stateInputEqConstraint.base.f.allFinite()) {
      errorDescription << "Input-state equality constraint is not finite.\n";
    }
    if (!data.stateInputEqConstraint.base.dfdx.allFinite()) {
      errorDescription << "Input-state equality constraint derivative w.r.t. state is not finite.\n";
    }
    if (!data.stateInputEqConstraint.base.dfdu.allFinite()) {
      errorDescription << "Input-state equality constraint derivative w.r.t. input is not finite.\n";
    }
    if (numConstraints > inputDim) {
      errorDescription << "Number of active state-input equality constraints (a.k.a. " + std::to_string(numConstraints) +
                              ") should be less-equal to the input dimension (a.k.a. " + std::to_string(inputDim) + ").\n";
    }
    const size_t DmRank = LinearAlgebra::rank(data.stateInputEqConstraint.base.dfdu);
    if (DmRank != numConstraints) {
      errorDescription << "Input-state equality constraint derivative w.r.t. input is not full-row rank. It's rank is " + std::to_string(DmRank) +
                              " while the expected rank is " + std::to_string(numConstraints) + ".\n";
    }
  }

  return errorDescription.str();
}

}  // namespace sto_ipm
}  // namespace ocs2
