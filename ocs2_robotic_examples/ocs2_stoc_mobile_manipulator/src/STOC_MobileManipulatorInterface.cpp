#include <string>

#include <pinocchio/fwd.hpp>  // forward declarations must be included first.
#include <pinocchio/multibody/joint/joint-composite.hpp>
#include <pinocchio/multibody/model.hpp>

#include <ocs2_core/misc/LoadData.h>
#include <ocs2_core/misc/LoadStdVectorOfPair.h>

#include "ocs2_stoc_mobile_manipulator/STOC_MobileManipulatorInterface.h"

// Boost
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

namespace ocs2 {
namespace mobile_manipulator {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
STOC_MobileManipulatorInterface::STOC_MobileManipulatorInterface(const std::string& taskFile, const std::string& libraryFolder,
                                                                 const std::string& urdfFile) 
 : baseInterface_(taskFile, libraryFolder, urdfFile) { 

  stocSettings_ = stoc::loadSettings(taskFile, "stoc");

  problem_ = ipm::OptimalControlProblem(baseInterface_.getOptimalControlProblem());
  const bool activateJointPositionLimit = problem_.softConstraintPtr->erase("jointPositionLimits");
  if (activateJointPositionLimit) {
    problem_.stateInequalityConstraintPtr->add("jointPositionLimits",
                                               getJointPositionLimitConstraint(getPinocchioInterface(), taskFile, "jointPositionLimits"));
  }
  problem_.softConstraintPtr->erase("jointVelocityLimits");
  problem_.inequalityConstraintPtr->add("jointVelocityLimits", getJointVelocityLimitConstraint(taskFile, "jointVelocityLimits"));
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::unique_ptr<StateConstraint> STOC_MobileManipulatorInterface::getJointPositionLimitConstraint(const PinocchioInterface& pinocchioInterface,
                                                                                                  const std::string& taskFile,
                                                                                                  const std::string& prefix) {
  const int armStateDim = getManipulatorModelInfo().armDim;
  const auto& model = pinocchioInterface.getModel();

  // arm joint DOF limits from the parsed URDF
  vector_t lowerBound = model.lowerPositionLimit.tail(armStateDim);
  vector_t upperBound = model.upperPositionLimit.tail(armStateDim);

  boost::property_tree::ptree pt;
  boost::property_tree::read_info(taskFile, pt);
  std::cerr << "\n #### JointPositionLimits Settings: ";
  std::cerr << "\n #### =============================================================================\n";
  std::cerr << " #### lowerBound: " << lowerBound.transpose() << '\n';
  std::cerr << " #### upperBound: " << upperBound.transpose() << '\n';
  std::cerr << " #### =============================================================================\n";

  std::unique_ptr<StateConstraint> constraint(new JointPositionLimitsIneq(armStateDim, lowerBound, upperBound));
  return constraint;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::unique_ptr<StateInputConstraint> STOC_MobileManipulatorInterface::getJointVelocityLimitConstraint(const std::string& taskFile,
                                                                                                       const std::string& prefix) {
  const auto& manipulatorModelInfo = getManipulatorModelInfo();
  const int baseInputDim = manipulatorModelInfo.inputDim - manipulatorModelInfo.armDim;
  const int armInputDim = manipulatorModelInfo.armDim;
  vector_t lowerBound = vector_t::Zero(manipulatorModelInfo.inputDim);
  vector_t upperBound = vector_t::Zero(manipulatorModelInfo.inputDim);

  // arm base DOFs velocity limits
  if (baseInputDim > 0) {
    vector_t lowerBoundBase = vector_t::Zero(baseInputDim);
    vector_t upperBoundBase = vector_t::Zero(baseInputDim);

    loadData::loadEigenMatrix(taskFile, prefix + ".lowerBound.base." + modelTypeEnumToString(manipulatorModelInfo.manipulatorModelType),
                              lowerBoundBase);
    loadData::loadEigenMatrix(taskFile, prefix + ".upperBound.base." + modelTypeEnumToString(manipulatorModelInfo.manipulatorModelType),
                              upperBoundBase);
    lowerBound.head(baseInputDim) = lowerBoundBase;
    upperBound.head(baseInputDim) = upperBoundBase;
  }

  // arm joint DOFs velocity limits
  vector_t lowerBoundArm = vector_t::Zero(armInputDim);
  vector_t upperBoundArm = vector_t::Zero(armInputDim);
  loadData::loadEigenMatrix(taskFile, prefix + ".lowerBound.arm", lowerBoundArm);
  loadData::loadEigenMatrix(taskFile, prefix + ".upperBound.arm", upperBoundArm);
  lowerBound.tail(armInputDim) = lowerBoundArm;
  upperBound.tail(armInputDim) = upperBoundArm;

  boost::property_tree::ptree pt;
  boost::property_tree::read_info(taskFile, pt);
  std::cerr << "\n #### JointVelocityLimits Settings: ";
  std::cerr << "\n #### =============================================================================\n";
  std::cerr << " #### 'lowerBound':  " << lowerBound.transpose() << std::endl;
  std::cerr << " #### 'upperBound':  " << upperBound.transpose() << std::endl;
  std::cerr << " #### =============================================================================\n";

  std::unique_ptr<StateInputConstraint> constraint(new JointVelocityLimitsIneq(manipulatorModelInfo.inputDim, 
                                                                               lowerBound, upperBound));
  return constraint;
}

}  // namespace mobile_manipulator
}  // namespace ocs2
