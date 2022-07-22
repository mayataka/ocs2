#include <iostream>
#include <string>

#include "ocs2_stoc_legged_robot/STOC_LeggedRobotInterface.h"

// Boost
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

namespace ocs2 {
namespace legged_robot {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
STOC_LeggedRobotInterface::STOC_LeggedRobotInterface(const std::string& taskFile, const std::string& urdfFile, const std::string& referenceFile)
  : baseInterface_(taskFile, urdfFile, referenceFile) {
  bool verbose;
  loadData::loadCppDataType(taskFile, "legged_robot_interface.verbose", verbose);
  stocSettings_ = stoc::loadSettings(taskFile, "stoc", verbose);
  setupOptimalConrolProblem(taskFile, urdfFile, referenceFile, verbose);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void STOC_LeggedRobotInterface::setupOptimalConrolProblem(const std::string& taskFile, const std::string& urdfFile,
                                                          const std::string& referenceFile, bool verbose) {
  problemPtr_ = std::unique_ptr<ipm::OptimalControlProblem>(new ipm::OptimalControlProblem(baseInterface_.getOptimalControlProblem()));
  const scalar_t frictionCoefficient = loadFrictionConeSettings(taskFile, verbose);
  const auto& centroidalModelInfo = baseInterface_.getCentroidalModelInfo();
  const auto& modelSettings = baseInterface_.modelSettings();
  for (size_t i = 0; i < centroidalModelInfo.numThreeDofContacts; i++) {
    const std::string& footName = modelSettings.contactNames3DoF[i];
    problemPtr_->softConstraintPtr->erase(footName + "_frictionCone");
    problemPtr_->inequalityConstraintPtr->add(footName + "_frictionCone", getFrictionConeConstraint(i, frictionCoefficient));
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t STOC_LeggedRobotInterface::loadFrictionConeSettings(const std::string& taskFile, bool verbose) const {
  boost::property_tree::ptree pt;
  boost::property_tree::read_info(taskFile, pt);
  const std::string prefix = "frictionConeConstraint.";

  scalar_t frictionCoefficient = 1.0;
  if (verbose) {
    std::cerr << "\n #### Friction Cone Settings: ";
    std::cerr << "\n #### =============================================================================\n";
  }
  loadData::loadPtreeValue(pt, frictionCoefficient, prefix + "frictionCoefficient", verbose);
  if (verbose) {
    std::cerr << " #### =============================================================================\n";
  }

  return frictionCoefficient;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::unique_ptr<StateInputConstraint> STOC_LeggedRobotInterface::getFrictionConeConstraint(size_t contactPointIndex, 
                                                                                           scalar_t frictionCoefficient) {
  FrictionConeIneqConstraint::Config frictionConeConConfig(frictionCoefficient);
  CentroidalModelInfo centroidalModelInfo = baseInterface_.getCentroidalModelInfo();
  std::unique_ptr<FrictionConeIneqConstraint> frictionConeConstraintPtr(
      new FrictionConeIneqConstraint(*baseInterface_.getSwitchedModelReferenceManagerPtr(), std::move(frictionConeConConfig), 
                                     contactPointIndex, std::move(centroidalModelInfo)));
  return frictionConeConstraintPtr;
}
 
}  // namespace legged_robot
}  // namespace ocs2
 