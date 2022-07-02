#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_oc/oc_data/PrimalSolution.h>
#include <ocs2_ipm/oc_data/DualVariable.h>
#include <ocs2_ipm/model_data/ModelData.h>

#include <ocs2_core/control/LinearController.h>

namespace ocs2 {
namespace stoc {

/**
 * Primal data container
 *
 * The design philosophy behind is to keep all member variables consistent. All (time, post, .., modelDataEventTime)
 * trajectories should be the rollout result of the controller
 *
 * There is one exception that breaks the consistency. When using an external controller to initialize the controller, it is obvious that
 * the rest of member variables are not the result of the controller. But they will be cleared and populated when runInit is called.
 */
struct PrimalDataContainer {
  PrimalSolution primalSolution;
  // intermediate model data trajectory
  std::vector<ipm::ModelData> modelDataTrajectory;
  // event times model data
  std::vector<ipm::ModelData> modelDataEventTimes;

  void swap(PrimalDataContainer& other) {
    primalSolution.swap(other.primalSolution);
    modelDataTrajectory.swap(other.modelDataTrajectory);
    modelDataEventTimes.swap(other.modelDataEventTimes);
  }

  void clear() {
    primalSolution.clear();
    modelDataTrajectory.clear();
    modelDataEventTimes.clear();
  }
};

/**
 * Dual data container
 *
 * The design philosophy behind is to keep all member variables consistent. valueFunctionTrajectory is the direct result of
 * (projectedModelData,riccatiModification) trajectories.
 *
 */
struct DualDataContainer {
  std::vector<ipm::DualVariable> dualVariableTrajectory;
  std::vector<ipm::DualVariableDirection> dualDirectionTrajectory;

  // Riccati solution coefficients
  std::vector<ScalarFunctionQuadraticApproximation> valueFunctionTrajectory;

  inline void swap(DualDataContainer& other) {
    dualVariableTrajectory.swap(other.dualVariableTrajectory);
    dualDirectionTrajectory.swap(other.dualDirectionTrajectory);
    valueFunctionTrajectory.swap(other.valueFunctionTrajectory);
  }

  inline void clear() {
    dualVariableTrajectory.clear();
    dualDirectionTrajectory.clear();
    valueFunctionTrajectory.clear();
  }
};

}  // namespace stoc
}  // namespace ocs2
