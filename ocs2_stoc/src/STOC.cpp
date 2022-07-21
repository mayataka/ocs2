#include "ocs2_stoc/STOC.h"

#include <ocs2_sqp/MultipleShootingInitialization.h>
#include <ocs2_core/control/FeedforwardController.h>
#include <ocs2_core/control/LinearController.h>

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
  for (size_t w = 0; w < settings.nThreads; w++) {
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
  primalData_.clear();
  ipmData_.clear();
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
    infoStream << "STOC Benchmarking\t   :\tAverage time [ms]   (% of total runtime)\n";
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

const std::vector<ipm::PerformanceIndex>& STOC::getIpmIterationsLog() const {
  if (ipmPerformanceIndeces_.empty()) {
    throw std::runtime_error("[STOC]: No performance log yet, no problem solved yet?");
  } else {
    return ipmPerformanceIndeces_;
  }
}

void STOC::runImpl(scalar_t initTime, const vector_t& initState, scalar_t finalTime) {
  if (settings_.printSolverStatus || settings_.printLinesearch) {
    std::cerr << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    std::cerr << "\n+++++++++++++ STOC solver is initialized ++++++++++++++";
    std::cerr << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  }

  // Initialize references
  for (auto& optimalControlProblem : optimalControlProblemStock_) {
    const auto& targetTrajectories = this->getReferenceManager().getTargetTrajectories();
    optimalControlProblem.targetTrajectoriesPtr = &targetTrajectories;
  }

  // Determine time discretization, taking into account event times.
  const auto& modeSchedule = this->getReferenceManager().getModeSchedule();
  const auto initTimeDiscretization = multiPhaseTimeDiscretizationGrid(initTime, finalTime, settings_.dt, modeSchedule);
  const auto numPhases = initTimeDiscretization.back().phase;

  std::cout << "initTimeDiscretization" << initTimeDiscretization << std::endl; 

  // Initialize the state, input, and Ipm variables
  vector_array_t stateTrajectory, inputTrajectory;
  initializeStateInputTrajectories(initTimeDiscretization, initState, stateTrajectory, inputTrajectory);
  vector_array_t costateTrajectory;
  initializeCostateTrajectories(initTimeDiscretization, stateTrajectory, inputTrajectory, costateTrajectory);
  std::vector<ipm::IpmVariables> ipmVariablesTrajectory;
  initializeIpmVariablesTrajectories(initTimeDiscretization, initState, stateTrajectory, inputTrajectory, ipmVariablesTrajectory);

  // Bookkeeping
  performanceIndeces_.clear();
  ipmPerformanceIndeces_.clear();

  // Directions
  vector_array_t dx, du, dlmd;
  scalar_array_t dts(numPhases, 0.0);
  std::vector<ipm::IpmVariablesDirection> ipmVariablesDirectionTrajectory;

  int iter = 0;
  auto convergence = Convergence::FALSE;
  while (convergence == Convergence::FALSE) {
    if (settings_.printSolverStatus || settings_.printLinesearch) {
      std::cerr << "\nSTOC iteration: " << iter << "\n";
    }
    // Make QP approximation of nonlinear primal-dual interior point method
    const auto timeDiscretization = multiPhaseTimeDiscretizationGrid(initTime, finalTime, settings_.dt, modeSchedule);
    linearQuadraticApproximationTimer_.startTimer();
    const auto performanceIndex = approximateOptimalControlProblem(timeDiscretization, initState, stateTrajectory, inputTrajectory, 
                                                                   costateTrajectory, ipmVariablesTrajectory);
    ipmPerformanceIndeces_.push_back(performanceIndex);
    performanceIndeces_.push_back(convert(performanceIndex));
    linearQuadraticApproximationTimer_.endTimer();

    // Reserve and resize directions
    const auto N = timeDiscretization.size() - 1;
    dx.reserve(N + 1); du.reserve(N); dlmd.reserve(N + 1); 
    while (dx.size() < N + 1) { 
      dx.push_back(vector_t::Zero(stateTrajectory[0].size())); 
    }
    while (du.size() < N) { 
      du.push_back(vector_t::Zero(inputTrajectory[0].size())); 
    }
    while (dlmd.size() < N + 1) { 
      dlmd.push_back(vector_t::Zero(costateTrajectory[0].size())); 
    }

    // Solve QP
    riccatiRecursionTimer_.startTimer();
    riccatiRecursion_.backwardRecursion(timeDiscretization, primalData_.modelDataTrajectory);
    dx[0] = initState - stateTrajectory[0];
    riccatiRecursion_.forwardRecursion(timeDiscretization, primalData_.modelDataTrajectory, dx, du, dlmd, dts);
    riccatiRecursionTimer_.endTimer();

    // Debugging 
    // for (int i=0; i<N; ++i) {
    //   std::cout << "dx[" << i << "]: " << dx[i].transpose() << "\n";
    //   std::cout << "du[" << i << "]: " << du[i].transpose() << "\n";
    //   std::cout << "dlmd[" << i << "]: " << dlmd[i].transpose() << "\n";
    // }
    // std::cout << "dx[" << N << "]: " << dx[N].transpose() << "\n";
    // std::cout << "dlmd[" << N << "]: " << dlmd[N].transpose() << "\n";

    // Select step sizes
    ipmVariablesDirectionTrajectory.reserve(N + 1);
    linesearchTimer_.startTimer();
    const auto stepSizes = selectStepSizes(timeDiscretization, ipmVariablesTrajectory, dx, du, dlmd, dts, ipmVariablesDirectionTrajectory);
    const scalar_t primalStepSize = stepSizes.primalStepSize;
    const scalar_t dualStepSize = stepSizes.dualStepSize;
    linesearchTimer_.endTimer();
    if (settings_.printLinesearch) {
      std::cerr << "[Linesearch] Primal step size: " << primalStepSize << ", Dual step size: " << dualStepSize << "\n";
    }

    // Update iterate
    updateIterate(stateTrajectory, inputTrajectory, costateTrajectory, ipmVariablesTrajectory, 
                  dx, du, dlmd, ipmVariablesDirectionTrajectory, primalStepSize, dualStepSize);

    // Check convergence
    convergence = checkConvergence(iter, performanceIndex, primalStepSize, dualStepSize);

    // Next iteration
    ++iter;
    ++totalNumIterations_;
  }

  computeControllerTimer_.startTimer();
  const auto timeDiscretization = multiPhaseTimeDiscretizationGrid(initTime, finalTime, settings_.dt, modeSchedule);
  setPrimalSolution(timeDiscretization, std::move(stateTrajectory), std::move(inputTrajectory));
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

void STOC::initializeStateInputTrajectories(const std::vector<Grid>& timeDiscretization, const vector_t& initState, 
                                            vector_array_t& stateTrajectory, vector_array_t& inputTrajectory) {
  const size_t N = static_cast<int>(timeDiscretization.size()) - 1;  // // size of the input trajectory
  stateTrajectory.clear();
  stateTrajectory.reserve(N + 1);
  inputTrajectory.clear();
  inputTrajectory.reserve(N);
  auto& primalSolution = primalData_.primalSolution;

  // Determine till when to use the previous solution
  scalar_t interpolateStateTill = timeDiscretization.front().time;
  scalar_t interpolateInputTill = timeDiscretization.front().time;
  if (primalSolution.timeTrajectory_.size() >= 2) {
    interpolateStateTill = primalSolution.timeTrajectory_.back();
    interpolateInputTill = primalSolution.timeTrajectory_[primalSolution.timeTrajectory_.size() - 2];
  }

  // Initial state
  const scalar_t initTime = getIntervalStart(timeDiscretization[0]);
  if (initTime < interpolateStateTill) {
    stateTrajectory.push_back(
        LinearInterpolation::interpolate(initTime, primalSolution.timeTrajectory_, primalSolution.stateTrajectory_));
  } else {
    stateTrajectory.push_back(initState);
  }

  for (int i = 0; i < N; i++) {
    if (timeDiscretization[i].event == Grid::Event::PreEvent) {
      // Event Node
      inputTrajectory.push_back(vector_t());  // no input at event node
      stateTrajectory.push_back(multiple_shooting::initializeEventNode(timeDiscretization[i].time, stateTrajectory.back()));
    } else {
      // Intermediate node
      const scalar_t time = getIntervalStart(timeDiscretization[i]);
      const scalar_t nextTime = getIntervalEnd(timeDiscretization[i + 1]);
      vector_t input, nextState;
      if (time > interpolateInputTill || nextTime > interpolateStateTill) {  // Using initializer
        std::tie(input, nextState) =
            multiple_shooting::initializeIntermediateNode(*initializerPtr_, time, nextTime, stateTrajectory.back());
      } else {  // interpolate previous solution
        std::tie(input, nextState) = multiple_shooting::initializeIntermediateNode(primalSolution, time, nextTime, stateTrajectory.back());
      }
      inputTrajectory.push_back(std::move(input));
      stateTrajectory.push_back(std::move(nextState));
    }
  }
}

void STOC::initializeCostateTrajectories(const std::vector<Grid>& timeDiscretization, const vector_array_t& stateTrajectory, 
                                         const vector_array_t& inputTrajectory, vector_array_t& costateTrajectory) {
  const size_t N = static_cast<int>(timeDiscretization.size()) - 1;  // // size of the input trajectory
  costateTrajectory.clear();
  costateTrajectory.reserve(N + 1);
  // Final costate
  auto& optimalControlProblem = optimalControlProblemStock_[0];
  const scalar_t tN = getIntervalStart(timeDiscretization[N]);
  ipm::ModelData modelData;
  ipm::approximateFinalLQ(optimalControlProblem, tN, stateTrajectory[N], modelData);
  const vector_t finalCostate = - modelData.cost.dfdx;
  for (int i = 0; i <= N; i++) {
    costateTrajectory.push_back(finalCostate);
  }
}

void STOC::initializeIpmVariablesTrajectories(const std::vector<Grid>& timeDiscretization, const vector_t& initState, 
                                              const vector_array_t& stateTrajectory, const vector_array_t& inputTrajectory, 
                                              std::vector<ipm::IpmVariables>& ipmVariablesTrajectory, scalar_t barrier) {
  const size_t N = static_cast<int>(timeDiscretization.size()) - 1;
  auto& modelDataTrajectory = primalData_.modelDataTrajectory;
  modelDataTrajectory.clear();
  modelDataTrajectory.reserve(N + 1);
  ipmVariablesTrajectory.clear();
  ipmVariablesTrajectory.resize(N + 1);

  std::atomic_int timeIndex{0};
  auto parallelTask = [&](int workerId) {
    // Get worker specific resources
    auto& optimalControlProblem = optimalControlProblemStock_[workerId];

    int i = timeIndex++;
    while (i < N) {
      if (timeDiscretization[i].event == Grid::Event::PreEvent) {
        // Event node
        const scalar_t ti = getIntervalStart(timeDiscretization[i]);
        ipm::approximatePreJumpLQ(optimalControlProblem, ti, stateTrajectory[i], modelDataTrajectory[i]);
        ipm::initIpmVariables(modelDataTrajectory[i], ipmVariablesTrajectory[i], barrier);
      } else {
        // Normal, intermediate node
        const scalar_t ti = getIntervalStart(timeDiscretization[i]);
        const scalar_t dt = getIntervalDuration(timeDiscretization[i], timeDiscretization[i + 1]);
        ipm::approximateIntermediateLQ(optimalControlProblem, ti, stateTrajectory[i], inputTrajectory[i], modelDataTrajectory[i]);
        ipm::initIpmVariables(modelDataTrajectory[i], ipmVariablesTrajectory[i], barrier);
      }

      i = timeIndex++;
    }

    if (i == N) {  // Only one worker will execute this
      const scalar_t tN = getIntervalStart(timeDiscretization[N]);
      ipm::approximateFinalLQ(optimalControlProblem, tN, stateTrajectory[N], modelDataTrajectory[N]);
      ipm::initIpmVariables(modelDataTrajectory[i], ipmVariablesTrajectory[i], barrier);
    }
  };
  runParallel(std::move(parallelTask));
}

ipm::PerformanceIndex STOC::approximateOptimalControlProblem(const std::vector<Grid>& timeDiscretization, const vector_t& initState, 
                                                             const vector_array_t& stateTrajectory, const vector_array_t& inputTrajectory,
                                                             const vector_array_t& costateTrajectory, 
                                                             const std::vector<ipm::IpmVariables>& ipmVariablesTrajectory) {
  // Problem horizon
  const int N = static_cast<int>(timeDiscretization.size()) - 1;
  // create alias
  auto& modelDataTrajectory = primalData_.modelDataTrajectory;
  auto& ipmDataTrajectory = ipmData_.ipmDataTrajectory;

  modelDataTrajectory.clear();
  modelDataTrajectory.resize(timeDiscretization.size());
  ipmDataTrajectory.clear();
  ipmDataTrajectory.resize(timeDiscretization.size());

  std::vector<ipm::PerformanceIndex> performance(settings_.nThreads, ipm::PerformanceIndex());

  std::atomic_int timeIndex{0};
  auto parallelTask = [&](int workerId) {
    // Get worker specific resources
    auto& optimalControlProblem = optimalControlProblemStock_[workerId];
    ipm::PerformanceIndex workerPerformance;  // Accumulate performance in local variable
    // const bool projection = settings_.projectStateInputEqualityConstraints;

    int i = timeIndex++;
    while (i < N) {
      if (timeDiscretization[i].event == Grid::Event::PreEvent) {
        // Event node
        const scalar_t ti = getIntervalStart(timeDiscretization[i]);
        ipm::approximatePreJumpLQ(optimalControlProblem, ti, stateTrajectory[i], modelDataTrajectory[i]);
        ipm::discretizePreJumpLQ(stateTrajectory[i], stateTrajectory[i+1], costateTrajectory[i], costateTrajectory[i+1], 
                                 modelDataTrajectory[i]);
        ipm::eliminateIpmVariablesPreJumpLQ(ipmVariablesTrajectory[i], modelDataTrajectory[i], ipmDataTrajectory[i]);
        workerPerformance += fromModelData(modelDataTrajectory[i]);
        workerPerformance += fromIpmData(ipmDataTrajectory[i]);
      } else {
        // Normal, intermediate node
        const scalar_t ti = getIntervalStart(timeDiscretization[i]);
        const scalar_t dt = getIntervalDuration(timeDiscretization[i], timeDiscretization[i + 1]);
        ipm::approximateIntermediateLQ(optimalControlProblem, ti, stateTrajectory[i], inputTrajectory[i], modelDataTrajectory[i]);
        ipm::discretizeIntermediateLQ(dt, stateTrajectory[i], stateTrajectory[i+1], costateTrajectory[i], costateTrajectory[i+1], 
                                      modelDataTrajectory[i]);
        ipm::eliminateIpmVariablesIntermediateLQ(ipmVariablesTrajectory[i], modelDataTrajectory[i], ipmDataTrajectory[i]);
        workerPerformance += fromModelData(modelDataTrajectory[i]);
        workerPerformance += fromIpmData(ipmDataTrajectory[i]);
      }

      i = timeIndex++;
    }

    if (i == N) {  // Only one worker will execute this
      const scalar_t tN = getIntervalStart(timeDiscretization[N]);
      ipm::approximateFinalLQ(optimalControlProblem, tN, stateTrajectory[N], modelDataTrajectory[N]);
      ipm::discretizeFinalLQ(costateTrajectory[N], modelDataTrajectory[N]);
      ipm::eliminateIpmVariablesFinalLQ(ipmVariablesTrajectory[N], modelDataTrajectory[N], ipmDataTrajectory[N]);
      workerPerformance += fromModelData(modelDataTrajectory[N]);
      workerPerformance += fromIpmData(ipmDataTrajectory[N]);
    }

    // Accumulate! Same worker might run multiple tasks
    performance[workerId] += workerPerformance;
  };
  runParallel(std::move(parallelTask));

  // Account for init state in performance
  performance.front().dynamicsViolationSSE += (initState - stateTrajectory.front()).squaredNorm();

  // Sum performance of the threads
  auto totalPerformance = std::accumulate(std::next(performance.begin()), performance.end(), performance.front());
  return totalPerformance;
}

STOC::StepSizes STOC::selectStepSizes(const std::vector<Grid>& timeDiscretization, 
                                      const std::vector<ipm::IpmVariables>& ipmVariablesTrajectory, const vector_array_t& dx, 
                                      const vector_array_t& du, const vector_array_t& dlmd, const scalar_array_t& dts,
                                      std::vector<ipm::IpmVariablesDirection>& ipmVariablesDirectionTrajectory) {
  // Problem horizon
  const int N = static_cast<int>(timeDiscretization.size()) - 1;
  // create alias
  const auto& modelDataTrajectory = primalData_.modelDataTrajectory;
  const auto& ipmDataTrajectory = ipmData_.ipmDataTrajectory;

  scalar_array_t primalStepSize(settings_.nThreads, 0.0);
  scalar_array_t dualStepSize(settings_.nThreads, 0.0);

  std::atomic_int timeIndex{0};
  auto parallelTask = [&](int workerId) {
    // Get worker specific resources
    scalar_t workerPrimalStepSize = 1.0;
    scalar_t workerDualStepSize = 1.0;

    int i = timeIndex++;
    while (i < N) {
      if (timeDiscretization[i].event == Grid::Event::PreEvent) {
        // Event node
        ipm::retrivePreJumpIpmVariablesDirection(modelDataTrajectory[i], ipmDataTrajectory[i], ipmVariablesTrajectory[i], dx[i], 
                                                 ipmVariablesDirectionTrajectory[i]);
        workerPrimalStepSize = std::min(ipm::preJumpPrimalStepSize(ipmDataTrajectory[i], ipmVariablesTrajectory[i], ipmVariablesDirectionTrajectory[i]), 
                                        workerPrimalStepSize);
        workerDualStepSize = std::min(ipm::preJumpDualStepSize(ipmDataTrajectory[i], ipmVariablesTrajectory[i], ipmVariablesDirectionTrajectory[i]), 
                                      workerDualStepSize);
      } else {
        // Normal, intermediate node
        ipm::retriveIntermediateIpmVariablesDirection(modelDataTrajectory[i], ipmDataTrajectory[i], ipmVariablesTrajectory[i], dx[i], du[i],
                                                      ipmVariablesDirectionTrajectory[i]);
        workerPrimalStepSize = std::min(ipm::intermediatePrimalStepSize(ipmDataTrajectory[i], ipmVariablesTrajectory[i], ipmVariablesDirectionTrajectory[i]), 
                                        workerPrimalStepSize);
        workerDualStepSize = std::min(ipm::intermediateDualStepSize(ipmDataTrajectory[i], ipmVariablesTrajectory[i], ipmVariablesDirectionTrajectory[i]), 
                                      workerDualStepSize);
      }
      i = timeIndex++;
    }

    if (i == N) {  // Only one worker will execute this
        ipm::retriveFinalIpmVariablesDirection(modelDataTrajectory[N], ipmDataTrajectory[N], ipmVariablesTrajectory[N], dx[N], 
                                               ipmVariablesDirectionTrajectory[N]);
        workerPrimalStepSize = std::min(ipm::finalPrimalStepSize(ipmDataTrajectory[i], ipmVariablesTrajectory[i], ipmVariablesDirectionTrajectory[i]), 
                                        workerPrimalStepSize);
        workerDualStepSize = std::min(ipm::finalDualStepSize(ipmDataTrajectory[i], ipmVariablesTrajectory[i], ipmVariablesDirectionTrajectory[i]), 
                                      workerDualStepSize);
    }

    // Accumulate! Same worker might run multiple tasks
    primalStepSize[workerId] = workerPrimalStepSize;
    dualStepSize[workerId] = workerDualStepSize;
  };
  runParallel(std::move(parallelTask));

  StepSizes stepSizes;
  stepSizes.primalStepSize = *std::min_element(primalStepSize.begin(), primalStepSize.end());
  stepSizes.dualStepSize = *std::min_element(dualStepSize.begin(), dualStepSize.end());
  return stepSizes;
}


void STOC::updateIterate(vector_array_t& x, vector_array_t& u, vector_array_t& lmd, std::vector<ipm::IpmVariables>& ipmVariablesTrajectory,
                         const vector_array_t& dx, const vector_array_t& du, const vector_array_t& dlmd,
                         const std::vector<ipm::IpmVariablesDirection>& ipmVariablesDirectionTrajectory,
                         scalar_t primalStepSize, scalar_t dualStepSize) {
  const size_t N = static_cast<size_t>(x.size()) - 1;
  for (size_t i=0; i<=N; ++i) {
    x[i].noalias() += primalStepSize * dx[i];
  }
  for (size_t i=0; i<N; ++i) {
    u[i].noalias() += primalStepSize * du[i];
  }
  for (size_t i=0; i<=N; ++i) {
    lmd[i].noalias() += primalStepSize * dlmd[i];
  }
  for (size_t i=0; i<=N; ++i) {
    ipm::updateIpmVariables(ipmVariablesTrajectory[i], ipmVariablesDirectionTrajectory[i], primalStepSize, dualStepSize);
  }
}

void STOC::setPrimalSolution(const std::vector<Grid>& timeDiscretization, vector_array_t&& stateTrajectory, vector_array_t&& inputTrajectory) {
  const size_t N = static_cast<size_t>(timeDiscretization.size()) - 1;
  // Clear old solution
  auto& primalSolution = primalData_.primalSolution;
  primalSolution.clear();

  // Correct for missing inputs at PreEvents
  for (size_t i = 0; i < N + 1; ++i) {
    if (timeDiscretization[i].event == Grid::Event::PreEvent && i > 0) {
      inputTrajectory[i] = inputTrajectory[i - 1];
    }
  }

  // Compute feedback, before x and u are moved to primal solution
  vector_array_t uff;
  matrix_array_t controllerGain;
  if (settings_.useFeedbackPolicy) {
    // see doc/LQR_full.pdf for detailed derivation for feedback terms
    uff = inputTrajectory;  // Copy and adapt in loop
    controllerGain.reserve(N + 1);
    const auto LQRPolicies = riccatiRecursion_.getLQRPolicies();
    for (int i = 0; (i + 1) < timeDiscretization.size(); i++) {
      const auto& KMatrix = LQRPolicies[i].K;
      if (timeDiscretization[i].event == Grid::Event::PreEvent && i > 0) {
        uff[i] = uff[i - 1];
        controllerGain.push_back(controllerGain.back());
      } else {
        // Linear controller has convention u = uff + K * x;
        // We computed u = u'(t) + K (x - x'(t));
        // >> uff = u'(t) - K x'(t)
        // if (constraintsProjection_[i].f.size() > 0) {
        //   controllerGain.push_back(std::move(constraintsProjection_[i].dfdx));  // Steal! Don't use after this.
        //   controllerGain.back().noalias() += constraintsProjection_[i].dfdu * KMatrix;
        // } else {
          controllerGain.push_back(KMatrix);
        // }
        uff[i].noalias() -= controllerGain.back() * stateTrajectory[i];
      }
    }
    // Copy last one to get correct length
    uff.push_back(uff.back());
    controllerGain.push_back(controllerGain.back());
  }

  // Construct nominal time, state and input trajectories
  primalSolution.stateTrajectory_ = std::move(stateTrajectory);
  inputTrajectory.push_back(inputTrajectory.back());  // Repeat last input to make equal length vectors
  primalSolution.inputTrajectory_ = std::move(inputTrajectory);
  primalSolution.timeTrajectory_.reserve(timeDiscretization.size());
  for (size_t i = 0; i < timeDiscretization.size(); i++) {
    primalSolution.timeTrajectory_.push_back(timeDiscretization[i].time);
    if (timeDiscretization[i].event == Grid::Event::PreEvent) {
      primalSolution.postEventIndices_.push_back(i + 1);
    }
  }
  primalSolution.modeSchedule_ = this->getReferenceManager().getModeSchedule();

  // Assign controller
  if (settings_.useFeedbackPolicy) {
    primalSolution.controllerPtr_.reset(new LinearController(primalSolution.timeTrajectory_, std::move(uff), std::move(controllerGain)));
  } else {
    primalSolution.controllerPtr_.reset(new FeedforwardController(primalSolution.timeTrajectory_, primalSolution.inputTrajectory_));
  }

}

STOC::Convergence STOC::checkConvergence(int iteration, const ipm::PerformanceIndex& performanceIndex,
                                         scalar_t primalStepSize, scalar_t dualStepSize) const {
  const auto primalFeas = performanceIndex.dynamicsViolationSSE 
                          + performanceIndex.equalityConstraintsSSE
                          + performanceIndex.inequalityConstraintsSSE;
  const auto dualFeas = performanceIndex.dualDynamicsViolationSSE
                          + performanceIndex.dualControlViolationSSE
                          + performanceIndex.complementarySlacknessSSE;
  if ((iteration + 1) >= settings_.numIteration) {
    // Falied to converge because the next iteration would exceed the specified number of iterations
    return Convergence::MAXITERATIONS;
  } else if (primalStepSize < settings_.minPrimalStepSize || dualStepSize < settings_.minDualStepSize) {
    // Failed to converge because step dual size is below the specified minimum
    return Convergence::STEPSIZE;
  } else if (primalFeas < settings_.primalFeasTol && dualFeas < settings_.dualFeasTol) {
    // Converged because the KKT error is below the specified tolerance
    return Convergence::SUCCESS;
  } else {
    // None of the above convergence criteria were met -> not converged.
    return Convergence::FALSE;
  }
}

}  // namespace ocs2
