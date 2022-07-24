#pragma once

// OCS2
#include <ocs2_core/Types.h>
#include <ocs2_legged_robot/LeggedRobotInterface.h>

#include <ocs2_ipm/oc_problem/OptimalControlProblem.h>
#include <ocs2_stoc/STOC_Settings.h>
#include <ocs2_stoc_legged_robot/constraint/FrictionConeIneqConstraint.h>

/**
 * LeggedRobotInterface class
 * General interface for mpc implementation on the legged robot model
 */
namespace ocs2 {
namespace legged_robot {

class STOC_LeggedRobotInterface {
 public:
  /**
   * Constructor
   *
   * @throw Invalid argument error if input task file or urdf file does not exist.
   *
   * @param [in] stocSettingFile: The absolute path to the configuration file for the STOC solver.
   * @param [in] taskFile: The absolute path to the configuration file for the MPC.
   * @param [in] urdfFile: The absolute path to the URDF file for the robot.
   * @param [in] referenceFile: The absolute path to the reference configuration file.
   */
  STOC_LeggedRobotInterface(const std::string& stocSettingFile, const std::string& taskFile, const std::string& urdfFile, 
                            const std::string& referenceFile);

  ~STOC_LeggedRobotInterface() = default;

  const ipm::OptimalControlProblem& getOptimalControlProblem() const { return *problemPtr_; }

  const ModelSettings& modelSettings() const { return baseInterface_.modelSettings(); }
  const stoc::Settings& stocSettings() const { return stocSettings_; }
  const mpc::Settings& mpcSettings() const { return baseInterface_.mpcSettings(); }

  const vector_t& getInitialState() const { return baseInterface_.getInitialState(); }
  PinocchioInterface& getPinocchioInterface() { return baseInterface_.getPinocchioInterface(); }
  const CentroidalModelInfo& getCentroidalModelInfo() const { return baseInterface_.getCentroidalModelInfo(); }
  std::shared_ptr<SwitchedModelReferenceManager> getSwitchedModelReferenceManagerPtr() const { 
    return baseInterface_.getSwitchedModelReferenceManagerPtr(); 
  }

  const LeggedRobotInitializer& getInitializer() const { return baseInterface_.getInitializer(); }
  std::shared_ptr<ReferenceManagerInterface> getReferenceManagerPtr() const { return baseInterface_.getReferenceManagerPtr(); }

 private:
  void setupOptimalConrolProblem(const std::string& taskFile, const std::string& urdfFile, const std::string& referenceFile, bool verbose);

  scalar_t loadFrictionConeSettings(const std::string& taskFile, bool verbose) const;
  std::unique_ptr<StateInputConstraint> getFrictionConeConstraint(size_t contactPointIndex, scalar_t frictionCoefficient);

  LeggedRobotInterface baseInterface_;
  std::unique_ptr<ipm::OptimalControlProblem> problemPtr_;
  stoc::Settings stocSettings_;
};

}  // namespace legged_robot
}  // namespace ocs2
