#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_oc/oc_data/PrimalSolution.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_ipm/oc_data/IpmVariables.h>
#include <ocs2_ipm/oc_data/IpmData.h>

#include <ocs2_core/control/LinearController.h>

namespace ocs2 {
namespace stoc {

/**
 * Primal data container
 */
struct PrimalDataContainer {
  PrimalSolution primalSolution;
  vector_array_t costateTrajectory;
  std::vector<ipm::ModelData> modelDataTrajectory;
  std::vector<ipm::ModelData> projectedModelDataTrajectory;
  std::vector<VectorFunctionLinearApproximation> constraintProjection;

  void swap(PrimalDataContainer& other) {
    primalSolution.swap(other.primalSolution);
    costateTrajectory.swap(other.costateTrajectory);
    modelDataTrajectory.swap(other.modelDataTrajectory);
    projectedModelDataTrajectory.swap(other.projectedModelDataTrajectory);
    constraintProjection.swap(other.constraintProjection);
  }

  void clear() {
    primalSolution.clear();
    costateTrajectory.clear();
    modelDataTrajectory.clear();
    projectedModelDataTrajectory.clear();
    constraintProjection.clear();
  }
};

/**
 * Ipm data container
 */
struct IpmDataContainer {
  std::vector<ipm::IpmVariables> ipmVariablesTrajectory;
  std::vector<ipm::IpmData> ipmDataTrajectory;

  inline void swap(IpmDataContainer& other) {
    ipmVariablesTrajectory.swap(other.ipmVariablesTrajectory);
    ipmDataTrajectory.swap(other.ipmDataTrajectory);
  }

  inline void clear() {
    ipmVariablesTrajectory.clear();
    ipmDataTrajectory.clear();
  }
};

}  // namespace stoc
}  // namespace ocs2
