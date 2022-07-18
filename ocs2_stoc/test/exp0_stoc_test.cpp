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
  settings.printSolverStatus = false;
  settings.printSolverStatistics = false;
  settings.printLinesearch = false;
  settings.nThreads = 4;

  std::vector<double> initEventTimes{1.0};
  std::vector<size_t> modeSequence{0, 1};
  auto referenceManagerPtr = getExp0ReferenceManager(initEventTimes, modeSequence);
  auto problem = createExp0Problem(referenceManagerPtr);
  auto ipmProblem = 

  const scalar_t startTime = 0.0;
  const scalar_t finalTime = 2.0;
  const vector_t initState = (vector_t(STATE_DIM) << 0.0, 2.0).finished();

  auto initializerPtr = std::unique_ptr<Initializer>(new DefaultInitializer(INPUT_DIM));

  STOC stoc(settings, problem, initializerPtr);

  /******************************************************************************************************/
  /******************************************************************************************************/
  /******************************************************************************************************/
  // STOC
  // stoc::STOC solver(&timeTriggeredRollout, &systemDerivative, &systemConstraint, &systemCostFunction, &operatingTrajectories,
  //                                 slqSettings, modeScheduleManagerPtr, nullptr, gddpSettings, nlpSettings);

  // // run ocs2 using LQ method for computing the derivatives
  // ocs2.gddpSettings().useLQForDerivatives_ = true;
  // ocs2.run(startTime, initState, finalTime, partitioningTimes, initEventTimes);
  // Eigen::VectorXd optimizedEventTimes_LQ;
  // ocs2.getParameters(optimizedEventTimes_LQ);
  // double optimizedCost_LQ;
  // ocs2.getCost(optimizedCost_LQ);

  // // run ocs2 using BVP method for computing the derivatives
  // ocs2.gddpSettings().useLQForDerivatives_ = false;
  // ocs2.run(startTime, initState, finalTime, partitioningTimes, initEventTimes);
  // Eigen::VectorXd optimizedEventTimes_BVP;
  // ocs2.getParameters(optimizedEventTimes_BVP);
  // double optimizedCost_BVP;
  // ocs2.getCost(optimizedCost_BVP);

  // /******************************************************************************************************/
  // /******************************************************************************************************/
  // /******************************************************************************************************/
  // std::cerr << "### Initial event times are:        ["
  //           << Eigen::Map<Eigen::VectorXd>(initEventTimes.data(), initEventTimes.size()).transpose() << "]\n";
  // std::cerr << "### Optimum cost LQ method:         " << optimizedCost_LQ << "\n";
  // std::cerr << "### Optimum event times LQ method:  [" << optimizedEventTimes_LQ.transpose() << "]\n";
  // std::cerr << "### Optimum cost BVP method:        " << optimizedCost_BVP << "\n";
  // std::cerr << "### Optimum event times BVP method: [" << optimizedEventTimes_BVP.transpose() << "]\n";

  // const double optimumCost = 9.766;
  // const std::vector<double> optimumEventTimes{0.1897};

  // ASSERT_NEAR(optimizedCost_LQ, optimumCost, 10 * slqSettings.ddpSettings_.minRelCost_)
  //     << "MESSAGE: OCS2 failed in the EXP1 using LQ approach for calculating derivatives!";
  // ASSERT_NEAR(optimizedCost_BVP, optimumCost, 10 * slqSettings.ddpSettings_.minRelCost_)
  //     << "MESSAGE: OCS2 failed in the EXP1 using BVP approach for calculating derivatives!";
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
