#include <gtest/gtest.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>

#include <ocs2_core/initialization/DefaultInitializer.h>
#include <ocs2_oc/test/circular_kinematics.h>

#include <ocs2_stoc/STOC.h>

using namespace ocs2;


TEST(test_circular_kinematics, solve_projected_EqConstraints) {
  // optimal control problem
  auto problem = createCircularKinematicsProblem("/tmp/sqp_test_generated");
  auto ipmProblem = ipm::OptimalControlProblem(problem);

  // Initializer
  DefaultInitializer zeroInitializer(2);

  // Solver settings
  stoc::Settings settings;
  settings.numIteration  = 10;
  settings.useFeedbackPolicy = true;
  settings.projectStateInputEqualityConstraints = true;
  settings.dt = 0.01;
  settings.printSolverStatus = true;
  settings.printSolverStatistics = true;
  settings.printLinesearch = true;
  settings.nThreads = 4;

  // Additional problem definitions
  const ocs2::scalar_t startTime = 0.0;
  const ocs2::scalar_t finalTime = 1.0;
  const ocs2::vector_t initState = (vector_t(2) << 1.0, 0.0).finished();  // radius 1.0

  // Solve
  STOC stoc(settings, ipmProblem, zeroInitializer);
  stoc.run(startTime, initState, finalTime);
  std::cout << stoc.getBenchmarkingInformation() << std::endl;
  std::cout << stoc.getIpmPerformanceIndeces() << std::endl;

  // Inspect solution
  const auto primalSolution = stoc.primalSolution(finalTime);
  for (int i = 0; i < primalSolution.timeTrajectory_.size(); i++) {
    std::cout << "time: " << std::setprecision(4) << primalSolution.timeTrajectory_[i] << "\t state: " << primalSolution.stateTrajectory_[i].transpose()
              << "\t input: " << primalSolution.inputTrajectory_[i].transpose() << std::endl;
  }

  // check constraint satisfaction
  for (int i = 0; i < primalSolution.timeTrajectory_.size()-1; i++) {
    if (primalSolution.inputTrajectory_[i].size() > 0) {
      EXPECT_NEAR(primalSolution.stateTrajectory_[i].dot(primalSolution.inputTrajectory_[i]), 0.0, 1.0e-06);
    }
  }

  // Check initial condition
  ASSERT_TRUE(primalSolution.stateTrajectory_.front().isApprox(initState));
  ASSERT_DOUBLE_EQ(primalSolution.timeTrajectory_.front(), startTime);
  ASSERT_DOUBLE_EQ(primalSolution.timeTrajectory_.back(), finalTime);

  // Check constraint satisfaction.
  const auto performance = stoc.getPerformanceIndeces();
  ASSERT_LT(performance.dynamicsViolationSSE, 1e-6);
  ASSERT_LT(performance.equalityConstraintsSSE, 1e-6);

  // Check feedback controller
  for (int i = 0; i < primalSolution.timeTrajectory_.size() - 1; i++) {
    const auto t = primalSolution.timeTrajectory_[i];
    const auto& x = primalSolution.stateTrajectory_[i];
    const auto& u = primalSolution.inputTrajectory_[i];
    // Feed forward part
    ASSERT_TRUE(u.isApprox(primalSolution.controllerPtr_->computeInput(t, x)));
  }
}