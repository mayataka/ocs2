#include <gtest/gtest.h>
#include <iostream>

#include <ocs2_centroidal_model/AccessHelperFunctions.h>
#include <ocs2_core/misc/LinearAlgebra.h>

#include "ocs2_stoc_legged_robot/constraint/FrictionConeIneqConstraint.h"
#include "ocs2_legged_robot/constraint/FrictionConeConstraint.h"
#include "ocs2_legged_robot/test/AnymalFactoryFunctions.h"

using namespace ocs2;
using namespace legged_robot;

class TestFrictionConeIneqConstraint : public testing::Test {
 public:
  using Matrix6x = Eigen::Matrix<scalar_t, 6, Eigen::Dynamic>;
  TestFrictionConeIneqConstraint() {}

  const PinocchioInterface pinocchioInterface = getPinocchioInterfaceFromUrdfFile(const std::string& urdfFile);
  const CentroidalModelType centroidalModelType = CentroidalModelType::SingleRigidBodyDynamics;
  std::unique_ptr<PinocchioInterface> pinocchioInterfacePtr = createAnymalPinocchioInterface();
  const CentroidalModelInfo centroidalModelInfo = createAnymalCentroidalModelInfo(pinocchioInterface, centroidalModelType);
  const auto nominalJointAngles
  const CentroidalModelInfo centroidalModelInfo createCentroidalModelInfo(pinocchioInterface, centroidalModelType, nominalJointAngles, 
                                                                          const std::vector<std::string>& threeDofContactNames, 
                                                                          {});

  const std::shared_ptr<SwitchedModelReferenceManager> referenceManagerPtr =
      createReferenceManager(centroidalModelInfo.numThreeDofContacts);
  PreComputation preComputation;
};


TEST_F(TestFrictionConeIneqConstraint, compareWithFrictionConeConstraint) {
  const FrictionConeConstraint::Config config;

  const size_t legNumber = 0;
  FrictionConeConstraint frictionConeConstraint(*referenceManagerPtr, config, legNumber, centroidalModelInfo);
  FrictionConeIneqConstraint frictionConeIneqConstraint(*referenceManagerPtr, config, legNumber, centroidalModelInfo);

  scalar_t t = 0.0;
  vector_t x = vector_t::Random(centroidalModelInfo.stateDim);
  vector_t u = vector_t::Random(centroidalModelInfo.inputDim);

  EXPECT_EQ(frictionConeConstraint.isActive(t), frictionConeIneqConstraint.isActive(t));
  EXPECT_EQ(frictionConeConstraint.getNumConstraints(t), frictionConeIneqConstraint.getNumConstraints(t));
  const auto linearApproximationSoftConstraint = -1 * frictionConeConstraint.getLinearApproximation(t, x, u, preComputation);
  const auto linearApproximation = frictionConeIneqConstraint.getLinearApproximation(t, x, u, preComputation);
  EXPECT_TRUE(linearApproximationSoftConstraint.f.isApprox(linearApproximation.f));
  EXPECT_TRUE(linearApproximationSoftConstraint.dfdx.isApprox(linearApproximation.dfdx));
  EXPECT_TRUE(linearApproximationSoftConstraint.dfdu.isApprox(linearApproximation.dfdu));
}
