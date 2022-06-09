#include "ocs2_stoc/STOC.h"

#include <iostream>
#include <numeric>

namespace ocs2 {

STOC::STOC(stoc::Settings settings, const ipm::OptimalControlProblem& optimalControlProblem, const Initializer& initializer)
    : SolverBase(),
      settings_(std::move(settings)),
      riccatiRecursion_(),
      threadPool_(std::max(settings_.nThreads, size_t(1)) - 1, settings_.threadPriority) {
  Eigen::setNbThreads(1);  // No multithreading within Eigen.
  Eigen::initParallel();

  // Clone objects to have one for each worker
  for (int w = 0; w < settings.nThreads; w++) {
    optimalControlProblemStock_.push_back(optimalControlProblem);
  }

  // Operating points
  initializerPtr_.reset(initializer.clone());

  // if (optimalControlProblem.equalityConstraintPtr->empty()) {
  //   settings_.projectStateInputEqualityConstraints = false;  // True does not make sense if there are no constraints.
  // }
}

STOC::~STOC() {
  if (settings_.printSolverStatistics) {
    std::cerr << getBenchmarkingInformation() << std::endl;
  }
}

void STOC::reset() {
  // Clear solution
  primalData_.primalSolution = PrimalSolution();
  performanceIndeces_.clear();

  // reset timers
  numProblems_ = 0;
  totalNumIterations_ = 0;
  linearQuadraticApproximationTimer_.reset();
  riccatiRecursionTimer_.reset();
  linesearchTimer_.reset();
  computeControllerTimer_.reset();
}

std::string STOC::getBenchmarkingInformation() const {
  const auto linearQuadraticApproximationTotal = linearQuadraticApproximationTimer_.getTotalInMilliseconds();
  const auto riccatiRecursionTotal = riccatiRecursionTimer_.getTotalInMilliseconds();
  const auto linesearchTotal = linesearchTimer_.getTotalInMilliseconds();
  const auto computeControllerTotal = computeControllerTimer_.getTotalInMilliseconds();

  const auto benchmarkTotal = linearQuadraticApproximationTotal + riccatiRecursionTotal + linesearchTotal + computeControllerTotal;

  std::stringstream infoStream;
  if (benchmarkTotal > 0.0) {
    const scalar_t inPercent = 100.0;
    infoStream << "\n########################################################################\n";
    infoStream << "The benchmarking is computed over " << totalNumIterations_ << " iterations. \n";
    infoStream << "SQP Benchmarking\t   :\tAverage time [ms]   (% of total runtime)\n";
    infoStream << "\tLQ Approximation   :\t" << linearQuadraticApproximationTimer_.getAverageInMilliseconds() << " [ms] \t\t("
               << linearQuadraticApproximationTotal / benchmarkTotal * inPercent << "%)\n";
    infoStream << "\tRiccati recursion  :\t" << riccatiRecursionTimer_.getAverageInMilliseconds() << " [ms] \t\t("
               << riccatiRecursionTotal / benchmarkTotal * inPercent << "%)\n";
    infoStream << "\tLinesearch         :\t" << linesearchTimer_.getAverageInMilliseconds() << " [ms] \t\t("
               << linesearchTotal / benchmarkTotal * inPercent << "%)\n";
    infoStream << "\tCompute Controller :\t" << computeControllerTimer_.getAverageInMilliseconds() << " [ms] \t\t("
               << computeControllerTotal / benchmarkTotal * inPercent << "%)\n";
  }
  return infoStream.str();
}

const std::vector<PerformanceIndex>& STOC::getIterationsLog() const {
  if (performanceIndeces_.empty()) {
    throw std::runtime_error("[STOC]: No performance log yet, no problem solved yet?");
  } else {
    return performanceIndeces_;
  }
}

void STOC::runImpl(scalar_t initTime, const vector_t& initState, scalar_t finalTime) {
  if (settings_.printSolverStatus || settings_.printLinesearch) {
    std::cerr << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    std::cerr << "\n+++++++++++++ STOC solver is initialized ++++++++++++++";
    std::cerr << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  }

  // Determine time discretization, taking into account event times.
  const auto& eventTimes = this->getReferenceManager().getModeSchedule().eventTimes;
  const auto timeDiscretization = stoc::multiPhaseTimeDiscretization(initTime, finalTime, settings_.dt, eventTimes);

  // // Initialize the state and input
  // vector_array_t x, u;
  // initializeStateInputTrajectories(initState, timeDiscretization, x, u);

  // Initialize references
  for (auto& optimalControlProblem : optimalControlProblemStock_) {
    const auto& targetTrajectories = this->getReferenceManager().getTargetTrajectories();
    optimalControlProblem.targetTrajectoriesPtr = &targetTrajectories;
  }

  // Bookkeeping
  performanceIndeces_.clear();

  int iter = 0;
  bool convergence = false;
  while (!convergence) {
    if (settings_.printSolverStatus || settings_.printLinesearch) {
      std::cerr << "\nSTOC iteration: " << iter << "\n";
    }
    // Make QP approximation
    linearQuadraticApproximationTimer_.startTimer();
    const auto performanceIndex = approximateOptimalControlProblem(initState, timeDiscretization);
    performanceIndeces_.push_back(performanceIndex);
    linearQuadraticApproximationTimer_.endTimer();

    // Solve QP
    riccatiRecursionTimer_.startTimer();
    const vector_t delta_x0 = initState - x[0];
    const auto deltaSolution = getOCPSolution(delta_x0);

    stoc::DiscreteTimeModeSchedule modeSchedule;
    riccatiRecursion_.backwardRecursion(modeSchedule, primalData_.modelDataTrajectory);
    vector_array_t dx;
    vector_array_t du;
    vector_array_t dlmd;
    scalar_array_t dts;
    riccatiRecursion_.forwardRecursion(modeSchedule, primalData_.modelDataTrajectory, dx, du, dlmd, dts);
    riccatiRecursionTimer_.endTimer();

    // Search step
    linesearchTimer_.startTimer();
    const auto stepSize = fractionToBoundaryRule(optimalControlProblemStock_, timeDiscretization, dx, du, dts, primalData_.modelDataTrajectory);
    linesearchTimer_.endTimer();

    // Check convergence
    convergence = checkConvergence(iter, baselinePerformance, stepInfo);

    // Next iteration
    ++iter;
    ++totalNumIterations_;
  }

  computeControllerTimer_.startTimer();
  setPrimalSolution(timeDiscretization, std::move(x), std::move(u));
  computeControllerTimer_.endTimer();

  ++numProblems_;

  if (settings_.printSolverStatus || settings_.printLinesearch) {
    std::cerr << "\nConvergence : " << toString(convergence) << "\n";
    std::cerr << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    std::cerr << "\n+++++++++++++ STOC solver has terminated ++++++++++++++";
    std::cerr << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  }
}

void STOC::runParallel(std::function<void(int)> taskFunction) {
  threadPool_.runParallel(std::move(taskFunction), settings_.nThreads);
}

PerformanceIndex STOC::approximateOptimalControlProblem(const vector_t& initState, 
                                                        const std::vector<AnnotatedTime>& timeDiscretization) {
  // Problem horizon
  const int N = static_cast<int>(timeDiscretization.size()) - 1;
  // create alias
  const auto& timeTrajectory = primalData_.primalSolution.timeTrajectory_;
  const auto& stateTrajectory = primalData_.primalSolution.stateTrajectory_;
  const auto& inputTrajectory = primalData_.primalSolution.inputTrajectory_;
  const auto& postEventIndices = primalData_.primalSolution.postEventIndices_;
  auto& modelDataTrajectory = primalData_.modelDataTrajectory;
  const auto& dualTrajectory = dualData_.dualVariableTrajectory;

  std::vector<PerformanceIndex> performance(settings_.nThreads, PerformanceIndex());
  // primalData_.resize() ???;

  std::atomic_int timeIndex{0};
  auto parallelTask = [&](int workerId) {
    // Get worker specific resources
    OptimalControlProblem& optimalControlProblem = optimalControlProblemStock_[workerId];
    PerformanceIndex workerPerformance;  // Accumulate performance in local variable
    // const bool projection = settings_.projectStateInputEqualityConstraints;

    int i = timeIndex++;
    while (i < N) {
      if (timeDiscretization[i].event == AnnotatedTime::Event::PreEvent) {
        // Event node
        const scalar_t ti = getIntervalStart(timeDiscretization[i]);
        ipm::approximatePreJumpLQ(optimalControlProblem, ti, stateTrajectory[i], inputTrajectory[i], modelDataTrajectory[i]);
        ipm::discretizePreJumpLQ(stateTrajectory[i], stateTrajectory[i+1],
                                 dualTrajectory[i], dualTrajectory[i+1], modelDataTrajectory[i]);
        workerPerformance += modelDataTrajectory[i].performanceIndex;
      } else {
        // Normal, intermediate node
        const scalar_t ti = getIntervalStart(timeDiscretization[i]);
        const scalar_t dt = getIntervalDuration(timeDiscretization[i], timeDiscretization[i + 1]);
        ipm::approximateIntermediateLQ(optimalControlProblem, ti, stateTrajectory[i], inputTrajectory[i], modelDataTrajectory[i]);
        ipm::discretizeIntermediateLQ(dt, stateTrajectory[i], stateTrajectory[i+1],
                                      dualTrajectory[i], dualTrajectory[i+1], modelDataTrajectory[i]);
        workerPerformance += modelDataTrajectory[i].performanceIndex;
      }

      i = timeIndex++;
    }

    if (i == N) {  // Only one worker will execute this
      const scalar_t tN = getIntervalStart(timeDiscretization[N]);
        ipm::approximateFinalLQ(optimalControlProblem, tN, stateTrajectory[N], modelDataTrajectory[N]);
        ipm::discretizeFinalLQ(dualTrajectory[N], modelDataTrajectory[N]);
      workerPerformance += modelDataTrajectory[N].performanceIndex;
    }

    // Accumulate! Same worker might run multiple tasks
    performance[workerId] += workerPerformance;
  };
  runParallel(std::move(parallelTask));

  // Account for init state in performance
  performance.front().dynamicsViolationSSE += (initState - x.front()).squaredNorm();

  // Sum performance of the threads
  PerformanceIndex totalPerformance = std::accumulate(std::next(performance.begin()), performance.end(), performance.front());
  totalPerformance.merit = totalPerformance.cost + totalPerformance.equalityLagrangian + totalPerformance.inequalityLagrangian;
  return totalPerformance;
}

}  // namespace ocs2
