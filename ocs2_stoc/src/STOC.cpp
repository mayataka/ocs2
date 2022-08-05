#include "ocs2_stoc/STOC.h"

#include <ocs2_sqp/MultipleShootingInitialization.h>
#include <ocs2_sqp/ConstraintProjection.h>
#include <ocs2_core/control/FeedforwardController.h>
#include <ocs2_core/control/LinearController.h>
#include <ocs2_core/misc/LinearInterpolation.h>
#include <ocs2_core/misc/Display.h>

#include <iostream>
#include <numeric>

namespace ocs2 {

STOC::STOC(stoc::Settings settings, const ipm::OptimalControlProblem& optimalControlProblem, const Initializer& initializer)
    : SolverBase(),
      settings_(std::move(settings)),
      riccatiRecursion_(settings_.riccatiSolverMode, settings_.switchingTimeTrustRegionRadius, settings_.enableSwitchingTimeTrustRegion),
      threadPool_(std::max(settings_.nThreads, size_t(1)) - 1, settings_.threadPriority) {
  Eigen::setNbThreads(1);  // No multithreading within Eigen.
  Eigen::initParallel();

  // Clone objects to have one for each worker
  for (size_t w = 0; w < settings.nThreads; w++) {
    optimalControlProblemStock_.push_back(optimalControlProblem);
  }

  // Operating points
  initializerPtr_.reset(initializer.clone());

  if (optimalControlProblem.equalityConstraintPtr->empty()) {
    settings_.projectStateInputEqualityConstraints = false;  // True does not make sense if there are no constraints.
  }

  if (optimalControlProblem.inequalityConstraintPtr->empty() 
      && optimalControlProblem.stateInequalityConstraintPtr->empty()
      && optimalControlProblem.preJumpInequalityConstraintPtr->empty()
      && optimalControlProblem.finalInequalityConstraintPtr->empty()) {
    settings_.targetBarrierParameter = settings_.initialBarrierParameter; // Turn off the barrier strategy if there are no inequality constraints.
  }
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
  ipmPerformanceIndeces_.clear();

  // reset timers
  numProblems_ = 0;
  totalNumIterations_ = 0;
  initializationTimer_.reset();
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

  // Save the initial mode schedule
  const auto initModeSchedule = this->getReferenceManager().getModeSchedule();

  // Determine time discretization, taking into account event times.
  auto modeSchedule = initModeSchedule;
  const auto initTimeDiscretization = multiPhaseTimeDiscretizationGrid(initTime, finalTime, settings_.dt, modeSchedule, 
                                                                       settings_.isStoEnabledInMode);
  const auto numPhases = initTimeDiscretization.back().phase + 1;

  // initialize Barrier param
  scalar_t barrierParameter = settings_.initialBarrierParameter;

  // Initialize the state, input, and costate
  vector_array_t stateTrajectory, inputTrajectory;
  initializeStateInputTrajectories(initTimeDiscretization, initState, stateTrajectory, inputTrajectory);
  vector_array_t costateTrajectory;
  initializeCostateTrajectories(initTimeDiscretization, stateTrajectory, inputTrajectory, costateTrajectory);

  // Ipm variables are initialized in approximateOptimalControlProblem() and approximateSwitchingTimeOptimizationProblem()
  std::vector<ipm::IpmVariables> ipmVariablesTrajectory;
  ipm::IpmVariables stoIpmVariables;

  // Bookkeeping
  performanceIndeces_.clear();
  ipmPerformanceIndeces_.clear();

  // Directions
  vector_array_t dx, du, dlmd;
  scalar_array_t dts(numPhases-1, 0.0);
  std::vector<ipm::IpmVariablesDirection> ipmVariablesDirectionTrajectory;
  ipm::IpmVariablesDirection stoIpmVariablesDirection;

  auto timeDiscretization = initTimeDiscretization;
  size_t iter = 0;
  bool initIpmVariables = true;
  auto convergence = Convergence::FALSE;
  while (convergence == Convergence::FALSE) {
    if (settings_.printSolverStatus || settings_.printLinesearch || settings_.printSwitchingTimeOptimization) {
      std::cerr << "\nSTOC iteration: " << iter << " (barrier parameter: " << barrierParameter << ")\n";
    }
    // Make QP approximation of nonlinear primal-dual interior point method
    linearQuadraticApproximationTimer_.startTimer();
    auto performanceIndex = approximateOptimalControlProblem(timeDiscretization, initState, stateTrajectory, inputTrajectory, 
                                                             costateTrajectory, ipmVariablesTrajectory, barrierParameter, 
                                                             initIpmVariables);
    performanceIndex += approximateSwitchingTimeOptimizationProblem(initTime, finalTime, initModeSchedule, modeSchedule,
                                                                    stoIpmVariables, barrierParameter, initIpmVariables);
    summarizeModelData(timeDiscretization);
    ipmPerformanceIndeces_.push_back(performanceIndex);
    performanceIndeces_.push_back(convert(performanceIndex));
    linearQuadraticApproximationTimer_.endTimer();

    // Reserve and resize directions
    const auto N = timeDiscretization.size() - 1;
    dx.resize(N + 1); du.resize(N); dlmd.resize(N + 1); 
    for (size_t i = 0; i < N + 1; ++i) {
      dx[i].resize(primalData_.modelDataTrajectory[i].stateDim);
    }
    for (size_t i = 0; i < N; ++i) {
      if (settings_.projectStateInputEqualityConstraints) {
        du[i].resize(primalData_.projectedModelDataTrajectory[i].inputDim);
      } else {
        du[i].resize(primalData_.modelDataTrajectory[i].inputDim);
      }
    }
    for (size_t i = 0; i < N + 1; ++i) {
      dlmd[i].resize(primalData_.modelDataTrajectory[i].stateDim);
    }

    // Solve QP
    riccatiRecursionTimer_.startTimer();
    if (settings_.projectStateInputEqualityConstraints) {
      riccatiRecursion_.backwardRecursion(timeDiscretization, primalData_.projectedModelDataTrajectory);
      dx[0] = initState - stateTrajectory[0];
      riccatiRecursion_.forwardRecursion(timeDiscretization, primalData_.projectedModelDataTrajectory, dx, du, dlmd, dts);
      vector_t tmp;  // 1 temporary for re-use.
      for (int i = 0; i < du.size(); i++) {
        const auto& constraintProjection = primalData_.constraintProjection[i];
        if (constraintProjection.f.size() > 0) {
          tmp.noalias() = constraintProjection.dfdu * du[i];
          du[i] = tmp + constraintProjection.f;
          du[i].noalias() += constraintProjection.dfdx * dx[i];
        }
      }
    } else {
      riccatiRecursion_.backwardRecursion(timeDiscretization, primalData_.modelDataTrajectory);
      dx[0] = initState - stateTrajectory[0];
      riccatiRecursion_.forwardRecursion(timeDiscretization, primalData_.modelDataTrajectory, dx, du, dlmd, dts);
    }
    riccatiRecursionTimer_.endTimer();

    // Select step sizes
    ipmVariablesDirectionTrajectory.reserve(N + 1);
    while (ipmVariablesDirectionTrajectory.size() < N + 1) {
      ipmVariablesDirectionTrajectory.push_back(ipm::IpmVariablesDirection());
    }
    linesearchTimer_.startTimer();
    const auto stepSizes = selectStepSizes(timeDiscretization, ipmVariablesTrajectory, dx, du, dlmd, ipmVariablesDirectionTrajectory,
                                           settings_.fractionToBoundaryMargin);
    const auto stepSizesSto = selectStepSizes(stoIpmVariables, dts, stoIpmVariablesDirection, settings_.fractionToBoundaryMargin);
    const scalar_t primalStepSize = std::min(stepSizes.primalStepSize, stepSizesSto.primalStepSize);
    const scalar_t dualStepSize = std::min(stepSizes.dualStepSize, stepSizesSto.dualStepSize);
    linesearchTimer_.endTimer();
    if (settings_.printLinesearch) {
      std::cerr << "[Linesearch] Primal step size: " << primalStepSize << ", Dual step size: " << dualStepSize << "\n";
    }

    // Update iterate
    updateIterate(stateTrajectory, inputTrajectory, costateTrajectory, ipmVariablesTrajectory, 
                  dx, du, dlmd, ipmVariablesDirectionTrajectory, primalStepSize, dualStepSize);
    updateIterate(initTime, finalTime, initModeSchedule, modeSchedule, stoIpmVariables, dts, stoIpmVariablesDirection, 
                  primalStepSize, dualStepSize);

    // Update mode schedule 
    // TODO: a safer way to update only the event times!
    getReferenceManager().setModeSchedule(modeSchedule);
    getReferenceManager().preSolverRun(initTime, finalTime, initState);

    updateTimeIntervals(initTime, finalTime, modeSchedule, timeDiscretization);
    const scalar_t maxTimeInterval = getMaxTimeInterval(timeDiscretization);
    if (settings_.printSwitchingTimeOptimization) {
      std::cerr << "[SwitchingTimeOptimization] Event times:   {" << toDelimitedString(modeSchedule.eventTimes) << "}\n";
    }

    // Check convergence
    convergence = checkConvergence(iter, barrierParameter, maxTimeInterval, performanceIndex, primalStepSize, dualStepSize);

    // mesh refinement
    if (maxTimeInterval > settings_.maxTimeInterval) {
      if (settings_.printSwitchingTimeOptimization) {
        std::cerr << "[SwitchingTimeOptimization] Mesh refinement is performed! Max time interval: " << maxTimeInterval << "\n";
      }
      setPrimalSolution(timeDiscretization, std::move(stateTrajectory), std::move(inputTrajectory), std::move(costateTrajectory));
      timeDiscretization = multiPhaseTimeDiscretizationGrid(initTime, finalTime, settings_.dt, modeSchedule, 
                                                            settings_.isStoEnabledInMode);
      initializeStateInputTrajectories(timeDiscretization, initState, stateTrajectory, inputTrajectory);
      initializeCostateTrajectories(timeDiscretization, stateTrajectory, inputTrajectory, costateTrajectory);
      initIpmVariables = true;
    }
    else {
      initIpmVariables = false;
    }

    // update barriere parameter
    barrierParameter = updateBarrierParameter(barrierParameter, performanceIndex);

    // Next iteration
    ++iter;
    ++totalNumIterations_;
  }

  computeControllerTimer_.startTimer();
  setPrimalSolution(timeDiscretization, std::move(stateTrajectory), std::move(inputTrajectory), std::move(costateTrajectory));
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
  const auto& primalSolution = primalData_.primalSolution;

  if (primalSolution.timeTrajectory_.size() >= 2) {
    // Linear interpolation if the previous solution is available
    for (size_t i = 0; i < N + 1; ++i) {
      costateTrajectory.push_back(
        LinearInterpolation::interpolate(timeDiscretization[i].time, primalSolution.timeTrajectory_, primalData_.costateTrajectory));
    }
  } else {
    // Fill with the final costate
    auto& optimalControlProblem = optimalControlProblemStock_[0];
    const scalar_t tN = getIntervalStart(timeDiscretization[N]);
    ipm::ModelData modelData;
    ipm::approximateFinalLQ(optimalControlProblem, tN, stateTrajectory[N], modelData);
    const vector_t finalCostate = - modelData.cost.dfdx;
    for (int i = 0; i <= N; i++) {
      costateTrajectory.push_back(finalCostate);
    }
  }
}

ipm::PerformanceIndex STOC::approximateOptimalControlProblem(const std::vector<Grid>& timeDiscretization, const vector_t& initState, 
                                                             const vector_array_t& stateTrajectory, const vector_array_t& inputTrajectory,
                                                             const vector_array_t& costateTrajectory, 
                                                             std::vector<ipm::IpmVariables>& ipmVariablesTrajectory, 
                                                             scalar_t barrierParameter, bool initIpmVariablesTrajectory) {
  // Problem horizon
  const int N = static_cast<int>(timeDiscretization.size()) - 1;
  // create alias
  auto& modelDataTrajectory = primalData_.modelDataTrajectory;
  auto& ipmDataTrajectory = ipmData_.ipmDataTrajectory;
  auto& projectedModelDataTrajectory = primalData_.projectedModelDataTrajectory;
  auto& constraintProjection = primalData_.constraintProjection;

  modelDataTrajectory.clear();
  modelDataTrajectory.resize(timeDiscretization.size());
  ipmDataTrajectory.clear();
  ipmDataTrajectory.resize(timeDiscretization.size());
  projectedModelDataTrajectory.clear();
  projectedModelDataTrajectory.resize(timeDiscretization.size());
  constraintProjection.clear();
  constraintProjection.resize(timeDiscretization.size());

  if (initIpmVariablesTrajectory) {
    ipmVariablesTrajectory.clear();
    ipmVariablesTrajectory.resize(N + 1);
  }

  const auto numGrids = getNumGrids(timeDiscretization);
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
        if (initIpmVariablesTrajectory) {
          ipm::initIpmVariables(modelDataTrajectory[i], ipmVariablesTrajectory[i], barrierParameter);
        }
        ipm::eliminateIpmVariablesPreJumpLQ(ipmVariablesTrajectory[i], modelDataTrajectory[i], ipmDataTrajectory[i], barrierParameter);
        if (settings_.projectStateInputEqualityConstraints) {
          projectedModelDataTrajectory[i] = modelDataTrajectory[i];
        }
        workerPerformance += ipm::fromModelData(modelDataTrajectory[i]);
        workerPerformance += ipm::fromIpmData(ipmDataTrajectory[i]);
      } else {
        // Normal, intermediate node
        const scalar_t ti = getIntervalStart(timeDiscretization[i]);
        const scalar_t dt = getIntervalDuration(timeDiscretization[i], timeDiscretization[i + 1]);
        const bool enableStateOnlyIneqConstraint = (i > 0); // State-only inequality constraints should be disabled at the initial time of the horizon.
        ipm::approximateIntermediateLQ(optimalControlProblem, ti, stateTrajectory[i], inputTrajectory[i], modelDataTrajectory[i],
                                       enableStateOnlyIneqConstraint);
        ipm::discretizeIntermediateLQ(dt, stateTrajectory[i], stateTrajectory[i+1], costateTrajectory[i], costateTrajectory[i+1], 
                                      modelDataTrajectory[i]);
        if (initIpmVariablesTrajectory) {
          ipm::initIpmVariables(modelDataTrajectory[i], ipmVariablesTrajectory[i], barrierParameter);
        }
        modelDataTrajectory[i].hamiltonian *= (1.0/static_cast<scalar_t>(numGrids[timeDiscretization[i].phase]));
        ipm::eliminateIpmVariablesIntermediateLQ(ipmVariablesTrajectory[i], modelDataTrajectory[i], ipmDataTrajectory[i], barrierParameter);
        if (settings_.projectStateInputEqualityConstraints) {
          ipm::projectIntermediateLQ(modelDataTrajectory[i], constraintProjection[i], projectedModelDataTrajectory[i]);
          workerPerformance += ipm::fromModelData(projectedModelDataTrajectory[i]);
        } else {
          workerPerformance += ipm::fromModelData(modelDataTrajectory[i]);
        }
        workerPerformance += ipm::fromIpmData(ipmDataTrajectory[i]);
      }

      i = timeIndex++;
    }

    if (i == N) {  // Only one worker will execute this
      const scalar_t tN = getIntervalStart(timeDiscretization[N]);
      ipm::approximateFinalLQ(optimalControlProblem, tN, stateTrajectory[N], modelDataTrajectory[N]);
      ipm::discretizeFinalLQ(costateTrajectory[N], modelDataTrajectory[N]);
      if (initIpmVariablesTrajectory) {
        ipm::initIpmVariables(modelDataTrajectory[N], ipmVariablesTrajectory[N], barrierParameter);
      }
      ipm::eliminateIpmVariablesFinalLQ(ipmVariablesTrajectory[N], modelDataTrajectory[N], ipmDataTrajectory[N], barrierParameter);
      if (settings_.projectStateInputEqualityConstraints) {
        projectedModelDataTrajectory[N] = modelDataTrajectory[N];
      }
      workerPerformance += ipm::fromModelData(modelDataTrajectory[N]);
      workerPerformance += ipm::fromIpmData(ipmDataTrajectory[N]);
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

ipm::PerformanceIndex STOC::approximateSwitchingTimeOptimizationProblem(scalar_t initTime, scalar_t finalTime, 
                                                                        const ModeSchedule& referenceModeSchedule, 
                                                                        const ModeSchedule& stoModeSchedule, 
                                                                        ipm::IpmVariables& stoIpmVariables,
                                                                        scalar_t barrierParameter, bool initStoIpmVariables) {
  // create alias
  auto& stoModelData = stoData_.stoModelData;
  auto& stoIpmData = ipmData_.stoIpmData;

  approximateStoProblem(optimalControlProblemStock_[0], initTime, finalTime, referenceModeSchedule, stoModeSchedule, stoModelData);
  if (initStoIpmVariables) {
    initIpmVariables(stoModelData, stoIpmVariables, barrierParameter);
  }
  eliminateIpmVariablesSTO(stoIpmVariables, stoModelData, stoIpmData, barrierParameter);

  const auto performance = ipm::fromStoModelData(stoModelData) + ipm::fromIpmData(stoIpmData);
  return performance;
}

void STOC::summarizeModelData(const std::vector<Grid>& timeDiscretization) {
  // Problem horizon
  const int N = static_cast<int>(timeDiscretization.size()) - 1;
  // create alias
  auto& modelDataTrajectory = primalData_.modelDataTrajectory;
  const auto& stoModelData = stoData_.stoModelData;

  for (size_t i = 0; i < N; ++i) {
    if (timeDiscretization[i].event == Grid::Event::PostEvent) {
      const auto switchIndex = timeDiscretization[i].phase - 1;
      modelDataTrajectory[i].hamiltonian.h += stoModelData.stoCost.dfdx.coeff(switchIndex);
      modelDataTrajectory[i].hamiltonian.dhdt = stoModelData.stoCost.dfdxx.coeff(switchIndex, switchIndex); // diagonal approximation
    }
  }
}

STOC::StepSizes STOC::selectStepSizes(const std::vector<Grid>& timeDiscretization, 
                                      const std::vector<ipm::IpmVariables>& ipmVariablesTrajectory, const vector_array_t& dx, 
                                      const vector_array_t& du, const vector_array_t& dlmd,
                                      std::vector<ipm::IpmVariablesDirection>& ipmVariablesDirectionTrajectory, 
                                      scalar_t fractionToBoundaryMargin) {
  // Problem horizon
  const int N = static_cast<int>(timeDiscretization.size()) - 1;
  // create alias
  const auto& modelDataTrajectory = primalData_.modelDataTrajectory;
  const auto& ipmDataTrajectory = ipmData_.ipmDataTrajectory;

  scalar_array_t primalStepSize(N+1, 1.0);
  scalar_array_t dualStepSize(N+1, 1.0);

  std::atomic_int timeIndex{0};
  auto parallelTask = [&](int workerId) {
    int i = timeIndex++;
    while (i < N) {
      if (timeDiscretization[i].event == Grid::Event::PreEvent) {
        // Event node
        ipm::retrivePreJumpIpmVariablesDirection(modelDataTrajectory[i], ipmDataTrajectory[i], ipmVariablesTrajectory[i], dx[i], 
                                                 ipmVariablesDirectionTrajectory[i]);
        primalStepSize[i] = ipm::preJumpPrimalStepSize(ipmDataTrajectory[i], ipmVariablesTrajectory[i], 
                                                       ipmVariablesDirectionTrajectory[i], fractionToBoundaryMargin);
        dualStepSize[i] = ipm::preJumpDualStepSize(ipmDataTrajectory[i], ipmVariablesTrajectory[i], 
                                                   ipmVariablesDirectionTrajectory[i], fractionToBoundaryMargin);
      } else {
        // Normal, intermediate node
        ipm::retriveIntermediateIpmVariablesDirection(modelDataTrajectory[i], ipmDataTrajectory[i], ipmVariablesTrajectory[i], dx[i], du[i],
                                                      ipmVariablesDirectionTrajectory[i]);
        primalStepSize[i] = ipm::intermediatePrimalStepSize(ipmDataTrajectory[i], ipmVariablesTrajectory[i], 
                                                            ipmVariablesDirectionTrajectory[i], fractionToBoundaryMargin);
        dualStepSize[i] = ipm::intermediateDualStepSize(ipmDataTrajectory[i], ipmVariablesTrajectory[i], 
                                                        ipmVariablesDirectionTrajectory[i], fractionToBoundaryMargin);
      }
      i = timeIndex++;
    }

    if (i == N) {  // Only one worker will execute this
        ipm::retriveFinalIpmVariablesDirection(modelDataTrajectory[N], ipmDataTrajectory[N], ipmVariablesTrajectory[N], dx[N], 
                                               ipmVariablesDirectionTrajectory[N]);
        primalStepSize[i] = ipm::finalPrimalStepSize(ipmDataTrajectory[N], ipmVariablesTrajectory[N],   
                                                     ipmVariablesDirectionTrajectory[N], fractionToBoundaryMargin);
        dualStepSize[i] = ipm::finalDualStepSize(ipmDataTrajectory[N], ipmVariablesTrajectory[N], 
                                                 ipmVariablesDirectionTrajectory[N], fractionToBoundaryMargin);
    }
  };
  runParallel(std::move(parallelTask));

  StepSizes stepSizes;
  stepSizes.primalStepSize = *std::min_element(primalStepSize.begin(), primalStepSize.end());
  stepSizes.dualStepSize = *std::min_element(dualStepSize.begin(), dualStepSize.end());
  return stepSizes;
}

STOC::StepSizes STOC::selectStepSizes(const ipm::IpmVariables& stoIpmVariables, const scalar_array_t& dts, 
                                      ipm::IpmVariablesDirection& stoIpmVariablesDirection, scalar_t fractionToBoundaryMargin) {
  // create alias
  const auto& stoModelData = stoData_.stoModelData;
  const auto& stoIpmData = ipmData_.stoIpmData;
  vector_t dts_(dts.size());
  for (size_t i=0; i<dts.size(); ++i) {
    dts_.coeffRef(i) = dts[i];
  }

  retriveStoIpmVariablesDirection(stoModelData, stoIpmData, stoIpmVariables, dts_, stoIpmVariablesDirection);
  StepSizes stepSizes;
  stepSizes.primalStepSize = stoPrimalStepSize(stoIpmData, stoIpmVariables, stoIpmVariablesDirection, fractionToBoundaryMargin);
  stepSizes.dualStepSize   = stoDualStepSize(stoIpmData, stoIpmVariables, stoIpmVariablesDirection, fractionToBoundaryMargin);
  return stepSizes;
}

void STOC::updateIterate(vector_array_t& x, vector_array_t& u, vector_array_t& lmd, std::vector<ipm::IpmVariables>& ipmVariablesTrajectory,
                         const vector_array_t& dx, const vector_array_t& du, const vector_array_t& dlmd,
                         const std::vector<ipm::IpmVariablesDirection>& ipmVariablesDirectionTrajectory,
                         scalar_t primalStepSize, scalar_t dualStepSize) {
  const size_t N = x.size() - 1;
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

void STOC::updateIterate(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, ModeSchedule& modeSchedule, 
                         ipm::IpmVariables& stoIpmVariables, const scalar_array_t& dts, 
                         const ipm::IpmVariablesDirection& stoIpmVariablesDirection, scalar_t primalStepSize, scalar_t dualStepSize) {
  if (settings_.isStoEnabledInMode.empty()) return;

  const auto validSwitchingTimeIndices = extractValidSwitchingTimeIndices(initTime, finalTime, referenceModeSchedule);
  for (size_t phase=0; phase<validSwitchingTimeIndices.size(); ++phase) {
    const auto mode = modeSchedule.modeSequence[phase+validSwitchingTimeIndices.front()];
    if (settings_.isStoEnabledInMode.find(mode) != settings_.isStoEnabledInMode.end()) {
      if (settings_.isStoEnabledInMode[mode]) modeSchedule.eventTimes[validSwitchingTimeIndices[phase]] += primalStepSize * dts[phase];
    }
  }
  ipm::updateIpmVariables(stoIpmVariables, stoIpmVariablesDirection, primalStepSize, dualStepSize);
}

void STOC::setPrimalSolution(const std::vector<Grid>& timeDiscretization, vector_array_t&& stateTrajectory, vector_array_t&& inputTrajectory,
                             vector_array_t&& costateTrajectory) {
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
        if (primalData_.constraintProjection[i].f.size() > 0) {
          controllerGain.push_back(std::move(primalData_.constraintProjection[i].dfdx));  // Steal! Don't use after this.
          controllerGain.back().noalias() += primalData_.constraintProjection[i].dfdu * KMatrix;
        } else {
          controllerGain.push_back(KMatrix);
        }
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

  // Construct costate trajectory
  primalData_.costateTrajectory = std::move(costateTrajectory);
}

STOC::Convergence STOC::checkConvergence(size_t iteration, scalar_t barrierParameter, scalar_t maxTimeInterval, 
                                         const ipm::PerformanceIndex& performanceIndex,
                                         scalar_t primalStepSize, scalar_t dualStepSize) const {
  const auto primalFeas = performanceIndex.dynamicsViolationSSE 
                          + performanceIndex.equalityConstraintsSSE
                          + performanceIndex.inequalityConstraintsSSE;
  const auto dualFeas = performanceIndex.dualDynamicsViolationSSE
                          + performanceIndex.dualControlViolationSSE
                          + performanceIndex.complementarySlacknessSSE;
  if (settings_.printSolverStatus) {
    std::cerr << "\nPrimal feasibility : " << primalFeas << ",  Dual feasibility : " << dualFeas << "\n";
  }
  if (primalFeas < settings_.primalFeasTol && dualFeas < settings_.dualFeasTol 
      && barrierParameter <= settings_.targetBarrierParameter && maxTimeInterval <= settings_.maxTimeInterval) {
    // Converged because the KKT error is below the specified tolerance and barrier parameter is below the target value
    return Convergence::SUCCESS;
  } else if (primalStepSize < settings_.minPrimalStepSize || dualStepSize < settings_.minDualStepSize) {
    // Failed to converge because step dual size is below the specified minimum
    return Convergence::STEPSIZE;
  } else if ((iteration + 1) >= settings_.numIteration) {
    // Falied to converge because the next iteration would exceed the specified number of iterations
    return Convergence::MAXITERATIONS;
  } else if (std::isnan(primalFeas) || std::isnan(dualFeas) || std::isinf(primalFeas) || std::isinf(dualFeas)) {
    // Falied to converge because the KKT error diverges
    return Convergence::DIVERGE;
  } else {
    // None of the above convergence criteria were met -> not converged.
    return Convergence::FALSE;
  }
}

scalar_t STOC::updateBarrierParameter(scalar_t currentBarrierParameter,
                                      const ipm::PerformanceIndex& performanceIndex) const {
  if (currentBarrierParameter <= settings_.targetBarrierParameter) {
    return currentBarrierParameter;
  }
  const auto primalFeas = performanceIndex.dynamicsViolationSSE 
                          + performanceIndex.equalityConstraintsSSE
                          + performanceIndex.inequalityConstraintsSSE;
  const auto dualFeas = performanceIndex.dualDynamicsViolationSSE
                          + performanceIndex.dualControlViolationSSE
                          + performanceIndex.complementarySlacknessSSE;
  if (primalFeas < settings_.barrierReductionPrimalFeasTol && dualFeas < settings_.barrierReductionDualFeasTol) {
    return std::min((currentBarrierParameter*settings_.barrierLinearDecreaseFactor), 
                    std::pow(currentBarrierParameter, settings_.barrierSuperlinearDecreasePower));
  }
  else {
    return currentBarrierParameter;
  }
}

}  // namespace ocs2
