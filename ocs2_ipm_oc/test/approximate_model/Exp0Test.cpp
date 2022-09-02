#include <gtest/gtest.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>

#include <ocs2_oc/test/EXP0.h>

#include <ocs2_ipm_oc/approximate_model/LinearQuadraticApproximator.h>
#include <ocs2_ipm_oc/approximate_model/LinearQuadraticDiscretizer.h>

using namespace ocs2;
using namespace ipm;

TEST(Exp0Test, Unconstrained) {
  static constexpr size_t STATE_DIM = 2;
  static constexpr size_t INPUT_DIM = 1;

  const scalar_array_t initEventTimes{1.0};
  const std::vector<size_t> modeSequence{0, 1};
  auto referenceManagerPtr = getExp0ReferenceManager(initEventTimes, modeSequence);
  auto problem = createExp0Problem(referenceManagerPtr);
  auto ipmProblem = ipm::OptimalControlProblem(problem);
  ipmProblem.targetTrajectoriesPtr = &referenceManagerPtr->getTargetTrajectories();

  const scalar_t time = 0.5;
  const scalar_t dt = 0.1;
  const vector_t state = vector_t::Random(STATE_DIM);
  const vector_t stateNext = vector_t::Random(STATE_DIM);
  const vector_t input = vector_t::Random(INPUT_DIM);
  const vector_t costate = vector_t::Random(STATE_DIM);
  const vector_t costateNext = vector_t::Random(STATE_DIM);

  const matrix_t A = (matrix_t(2, 2) << 0.6, 1.2, -0.8, 3.4).finished();
  const matrix_t B = (matrix_t(2, 1) << 1, 1).finished();

  ModelData intermediateModelData;
  approximateIntermediateLQ(ipmProblem, time, state, input, intermediateModelData);
  discretizeIntermediateLQ(dt, state, stateNext, costate, costateNext, intermediateModelData);
  const vector_t lx = dt * (matrix_t(2, 2) << 0.0, 0.0, 0.0, 1.0).finished() * (state - referenceManagerPtr->getTargetTrajectories().getDesiredState(time))
                      + (matrix_t::Identity(2, 2) + dt * A).transpose() * costateNext - costate;
  const vector_t lu = dt * (matrix_t::Identity(1, 1) * (input - referenceManagerPtr->getTargetTrajectories().getDesiredInput(time)))
                      + (dt * B).transpose() * costateNext;
  EXPECT_TRUE(intermediateModelData.cost.dfdxx.isApprox(dt * (matrix_t(2, 2) << 0.0, 0.0, 0.0, 1.0).finished()));
  EXPECT_TRUE(intermediateModelData.cost.dfduu.isApprox(dt* matrix_t::Identity(1, 1)));
  EXPECT_TRUE(intermediateModelData.cost.dfdux.isZero());
  EXPECT_TRUE(intermediateModelData.cost.dfdx.isApprox(lx));
  EXPECT_TRUE(intermediateModelData.cost.dfdu.isApprox(lu));
  EXPECT_TRUE(intermediateModelData.dynamics.f.isApprox(state+dt*(A*state+B*input)-stateNext));
  EXPECT_TRUE(intermediateModelData.dynamics.dfdx.isApprox(matrix_t::Identity(2, 2) + dt*A));
  EXPECT_TRUE(intermediateModelData.dynamics.dfdu.isApprox(dt*B));

  ModelData preJumpModelData;
  approximatePreJumpLQ(ipmProblem, initEventTimes[0], state, preJumpModelData);
  discretizePreJumpLQ(state, stateNext, costate, costateNext, preJumpModelData);
  EXPECT_TRUE(preJumpModelData.cost.dfdxx.isZero());
  EXPECT_TRUE(preJumpModelData.cost.dfduu.isZero());
  EXPECT_TRUE(preJumpModelData.cost.dfdux.isZero());
  EXPECT_TRUE(preJumpModelData.cost.dfdx.isApprox(costateNext-costate));
  EXPECT_TRUE(preJumpModelData.cost.dfdu.isZero());
  EXPECT_TRUE(preJumpModelData.dynamics.f.isApprox(state-stateNext));
  EXPECT_TRUE(preJumpModelData.dynamics.dfdx.isApprox(matrix_t::Identity(2, 2)));
  EXPECT_TRUE(preJumpModelData.dynamics.dfdu.isZero());

  ModelData finalModelData;
  const vector_t lxFi = matrix_t::Identity(2, 2) * (state - referenceManagerPtr->getTargetTrajectories().getDesiredState(time)) - costate;
  approximateFinalLQ(ipmProblem, time, state, finalModelData);
  discretizeFinalLQ(costate, finalModelData);
  EXPECT_TRUE(finalModelData.cost.dfdxx.isApprox(matrix_t::Identity(2, 2)));
  EXPECT_TRUE(finalModelData.cost.dfduu.isZero());
  EXPECT_TRUE(finalModelData.cost.dfdux.isZero());
  EXPECT_TRUE(finalModelData.cost.dfdx.isApprox(lxFi));
  EXPECT_TRUE(finalModelData.cost.dfdu.isZero());
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
