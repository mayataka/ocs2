#include <gtest/gtest.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>

#include <ocs2_core/initialization/DefaultInitializer.h>
#include <ocs2_oc/test/EXP0.h>

#include <ocs2_stoc/STOC.h>

using namespace ocs2;

TEST(exp0_stoc_test, exp0_stoc_test) {
  static constexpr size_t STATE_DIM = 2;
  static constexpr size_t INPUT_DIM = 1;

  stoc::Settings settings;
  settings.numIteration  = 10;
  settings.primalFeasTol = 1e-6;
  settings.dualFeasTol   = 1e-4;
  settings.initialBarrierParameter = 1.0e-01;
  settings.targetBarrierParameter = 1.0e-03;
  settings.barrierReductionRate = 0.5;
  settings.fractionToBoundaryMargin = 0.995;
  settings.useFeedbackPolicy = true;
  settings.dt = 0.01;
  settings.printSolverStatus = true;
  settings.printSolverStatistics = true;
  settings.printLinesearch = true;
  settings.nThreads = 1;

  std::vector<double> initEventTimes{1.0};
  std::vector<size_t> modeSequence{0, 1};
  auto referenceManagerPtr = getExp0ReferenceManager(initEventTimes, modeSequence);
  auto problem = createExp0Problem(referenceManagerPtr);
  auto ipmProblem = ipm::OptimalControlProblem(problem);

  const scalar_t startTime = 0.0;
  const scalar_t finalTime = 2.0;
  const vector_t initState = (vector_t(STATE_DIM) << 0.0, 2.0).finished();

  auto initializerPtr = std::unique_ptr<Initializer>(new DefaultInitializer(INPUT_DIM));

  STOC stoc(settings, ipmProblem, *initializerPtr);
  stoc.setReferenceManager(referenceManagerPtr);
  stoc.run(startTime, initState, finalTime);
  std::cout << stoc.getBenchmarkingInformation() << std::endl;
  std::cout << stoc.getPerformanceIndeces() << std::endl;
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
