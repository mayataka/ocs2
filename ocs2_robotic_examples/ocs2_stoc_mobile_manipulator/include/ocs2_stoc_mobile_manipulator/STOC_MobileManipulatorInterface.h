#pragma once

// OCS2
#include <ocs2_core/Types.h>
#include <ocs2_mobile_manipulator/MobileManipulatorInterface.h>

#include <ocs2_ipm/oc_problem/OptimalControlProblem.h>
#include <ocs2_stoc/STOC_Settings.h>
#include <ocs2_stoc_mobile_manipulator/constraint/JointPositionLimitsIneq.h>
#include <ocs2_stoc_mobile_manipulator/constraint/JointVelocityLimitsIneq.h>

namespace ocs2 {
namespace mobile_manipulator {

/**
 * Mobile Manipulator Robot Interface class
 */
class STOC_MobileManipulatorInterface {
 public:
  /**
   * Constructor
   *
   * @note Creates directory for generated library into if it does not exist.
   * @throw Invalid argument error if input task file or urdf file does not exist.
   *
   * @param [in] stocFile: The absolute path to the configuration file for the STOC solver.
   * @param [in] taskFile: The absolute path to the configuration file for the MPC.
   * @param [in] libraryFolder: The absolute path to the directory to generate CppAD library into.
   * @param [in] urdfFile: The absolute path to the URDF file for the robot.
   */
  STOC_MobileManipulatorInterface(const std::string& stocFile, const std::string& taskFile, const std::string& libraryFolder, 
                                  const std::string& urdfFile);

  const vector_t& getInitialState() { return baseInterface_.getInitialState(); }

  mpc::Settings& mpcSettings() { return baseInterface_.mpcSettings(); }

  stoc::Settings& stocSettings() { return stocSettings_; }

  const ipm::OptimalControlProblem& getOptimalControlProblem() const { return problem_; }

  std::shared_ptr<ReferenceManagerInterface> getReferenceManagerPtr() const { return baseInterface_.getReferenceManagerPtr(); }

  const Initializer& getInitializer() const { return baseInterface_.getInitializer(); }

  const RolloutBase& getRollout() const { return baseInterface_.getRollout(); }

  const PinocchioInterface& getPinocchioInterface() const { return baseInterface_.getPinocchioInterface(); }

  const ManipulatorModelInfo& getManipulatorModelInfo() const { return baseInterface_.getManipulatorModelInfo(); }

 private:
  std::unique_ptr<StateConstraint> getJointPositionLimitConstraint(const PinocchioInterface& pinocchioInterface, const std::string& taskFile,
                                                                   const std::string& prefix);
  std::unique_ptr<StateInputConstraint> getJointVelocityLimitConstraint(const std::string& taskFile, const std::string& prefix);


  MobileManipulatorInterface baseInterface_;

  stoc::Settings stocSettings_;

  ipm::OptimalControlProblem problem_;
};

}  // namespace mobile_manipulator
}  // namespace ocs2
