#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_core/control/LinearController.h>
#include <ocs2_core/initialization/Initializer.h>
#include <ocs2_core/misc/Benchmark.h>
#include <ocs2_core/misc/LinearInterpolation.h>
#include <ocs2_core/misc/Numerics.h>
#include <ocs2_core/thread_support/ThreadPool.h>

#include <ocs2_oc/oc_solver/SolverBase.h>

#include <ocs2_ipm/oc_problem/OptimalControlProblem.h>
#include <ocs2_ipm/oc_solver/PerformanceIndex.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_ipm/approximate_model/LinearQuadraticApproximator.h>
#include <ocs2_ipm/approximate_model/LinearQuadraticDiscretizer.h>
#include <ocs2_ipm/approximate_model/IpmVariablesEliminator.h>
#include <ocs2_ipm/approximate_model/IpmVariablesRetriver.h>
#include <ocs2_ipm/approximate_model/ConstraintProjection.h>

#include "ocs2_stoc/sto/LinearQuadraticApproximator.h"
#include "ocs2_stoc/sto/IpmVariables.h"
#include "ocs2_stoc/sto/IpmVariablesEliminator.h"
#include "ocs2_stoc/sto/IpmVariablesRetriver.h"

#include "ocs2_stoc/riccati_recursion/RiccatiRecursion.h"
#include "ocs2_stoc/TimeDiscretization.h"
#include "ocs2_stoc/STOC_Settings.h"
#include "ocs2_stoc/STOC_Data.h"

namespace ocs2 {

class STOC : public SolverBase {
 public:
  /**
   * Constructor
   *
   * @param settings : settings for the STOC solver.
   * @param [in] optimalControlProblem: The optimal control problem formulation.
   * @param [in] initializer: This class initializes the state-input for the time steps that no controller is available.
   */
  STOC(stoc::Settings settings, const ipm::OptimalControlProblem& optimalControlProblem, const Initializer& initializer);

  ~STOC() override;

  void reset() override;

  const PerformanceIndex& getPerformanceIndeces() const override { return getIterationsLog().back(); };

  size_t getNumIterations() const override { return totalNumIterations_; }

  const std::vector<PerformanceIndex>& getIterationsLog() const override;

  const ipm::PerformanceIndex& getIpmPerformanceIndeces() const { return getIpmIterationsLog().back(); };

  const std::vector<ipm::PerformanceIndex>& getIpmIterationsLog() const;

  scalar_t getFinalTime() const override { return primalData_.primalSolution.timeTrajectory_.back(); };

  void getPrimalSolution(scalar_t finalTime, PrimalSolution* primalSolutionPtr) const override { *primalSolutionPtr = primalData_.primalSolution; }

  ScalarFunctionQuadraticApproximation getValueFunction(scalar_t time, const vector_t& state) const override {
    throw std::runtime_error("[STOC] getValueFunction() not available yet.");
  };

  ScalarFunctionQuadraticApproximation getHamiltonian(scalar_t time, const vector_t& state, const vector_t& input) override {
    throw std::runtime_error("[STOC] getHamiltonian() not available yet.");
  }

  vector_t getStateInputEqualityConstraintLagrangian(scalar_t time, const vector_t& state) const override {
    throw std::runtime_error("[STOC] getStateInputEqualityConstraintLagrangian() not available yet.");
  }

  std::string getBenchmarkingInformation() const;

  /**
   * Const access to ddp settings
   */
  const stoc::Settings& settings() const { return settings_; }

 private:
  void runImpl(scalar_t initTime, const vector_t& initState, scalar_t finalTime) override;

  void runImpl(scalar_t initTime, const vector_t& initState, scalar_t finalTime, const ControllerBase* externalControllerPtr) override {
    if (externalControllerPtr == nullptr) {
      runImpl(initTime, initState, finalTime);
    } else {
      throw std::runtime_error("[STOC::run] This solver does not support external controller!");
    }
  }

  /** Run a task in parallel with settings.nThreads */
  void runParallel(std::function<void(int)> taskFunction);

  void initializeStateInputTrajectories(const std::vector<Grid>& timeDiscretization, const vector_t& initState,  
                                        vector_array_t& stateTrajectory, vector_array_t& inputTrajectory);

  void initializeCostateTrajectories(const std::vector<Grid>& timeDiscretization, const vector_array_t& stateTrajectory, 
                                     const vector_array_t& inputTrajectory, vector_array_t& costateTrajectory);

  ipm::PerformanceIndex approximateOptimalControlProblem(const std::vector<Grid>& timeDiscretization, const vector_t& initState, 
                                                         const vector_array_t& stateTrajectory, const vector_array_t& inputTrajectory,
                                                         const vector_array_t& costateTrajectory, 
                                                         std::vector<ipm::IpmVariables>& ipmVariablesTrajectory,
                                                         scalar_t barrierParameter, bool initIpmVariablesTrajectory);

  ipm::PerformanceIndex approximateSwitchingTimeOptimizationProblem(scalar_t initTime, scalar_t finalTime, 
                                                                    const ModeSchedule& referenceModeSchedule, 
                                                                    const ModeSchedule& stoModeSchedule, 
                                                                    ipm::IpmVariables& stoIpmVariables,
                                                                    scalar_t barrierParameter, bool initStoIpmVariables);

  void summarizeModelData(const std::vector<Grid>& timeDiscretization);

  struct StepSizes {
    scalar_t primalStepSize, dualStepSize;
  };
  StepSizes selectStepSizes(const std::vector<Grid>& timeDiscretization, const std::vector<ipm::IpmVariables>& ipmVariablesTrajectory,
                            const vector_array_t& dx, const vector_array_t& du, const vector_array_t& dlmd,
                            std::vector<ipm::IpmVariablesDirection>& ipmVariablesDirectionTrajectory, scalar_t fractionToBoundaryMargin); 

  StepSizes selectStepSizes(const ipm::IpmVariables& stoIpmVariables, const scalar_array_t& dts, 
                            ipm::IpmVariablesDirection& stoIpmVariablesDirection, scalar_t fractionToBoundaryMargin); 

  static void updateIterate(vector_array_t& x, vector_array_t& u, vector_array_t& lmd, std::vector<ipm::IpmVariables>& ipmVariablesTrajectory,
                            const vector_array_t& dx, const vector_array_t& du, const vector_array_t& dlmd, 
                            const std::vector<ipm::IpmVariablesDirection>& ipmVariablesDirectionTrajectory,
                            scalar_t primalStepSize, scalar_t dualStepSize);

  void updateIterate(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, ModeSchedule& modeSchedule, 
                     ipm::IpmVariables& stoIpmVariables, const scalar_array_t& dts, 
                     const ipm::IpmVariablesDirection& stoIpmVariablesDirection, scalar_t primalStepSize, scalar_t dualStepSize);

  void setPrimalSolution(const std::vector<Grid>& timeDiscretization, vector_array_t&& stateTrajectory, vector_array_t&& inputTrajectory,
                         vector_array_t&& costateTrajectory);

  enum class Convergence { FALSE, SUCCESS, MAXITERATIONS, STEPSIZE, DIVERGE };
  static std::string toString(Convergence convergence) {
    switch (convergence) {
    case Convergence::FALSE:
      return "FALSE";
      break;
    case Convergence::SUCCESS:
      return "SUCCESS";
      break;
    case Convergence::MAXITERATIONS:
      return "MAXITERATIONS";
      break;
    case Convergence::STEPSIZE:
      return "STEPSIZE";
      break;
    case Convergence::DIVERGE:
      return "DIVERGE";
      break;
    default:
      return "";
      break;
    }
  }
  Convergence checkConvergence(size_t iteration, scalar_t barrierParameter, const ipm::PerformanceIndex& performanceIndex, 
                               scalar_t primalStepSize, scalar_t dualStepSize) const;

  scalar_t updateBarrierParameter(scalar_t currentBarrierParameter,
                                  const ipm::PerformanceIndex& performanceIndex) const;
  
  // Problem definition
  stoc::Settings settings_;
  std::vector<ipm::OptimalControlProblem> optimalControlProblemStock_;
  std::unique_ptr<Initializer> initializerPtr_;

  // Data
  stoc::PrimalDataContainer primalData_;
  stoc::IpmDataContainer ipmData_;
  stoc::StoDataContainer stoData_;

  // riccati Recurision
  stoc::RiccatiRecursion riccatiRecursion_;

  // Threading
  ThreadPool threadPool_;

  // Iteration performance log
  std::vector<PerformanceIndex> performanceIndeces_;
  std::vector<ipm::PerformanceIndex> ipmPerformanceIndeces_;

  // Benchmarking
  size_t numProblems_{0};
  size_t totalNumIterations_{0};
  benchmark::RepeatedTimer initializationTimer_;
  benchmark::RepeatedTimer linearQuadraticApproximationTimer_;
  benchmark::RepeatedTimer riccatiRecursionTimer_;
  benchmark::RepeatedTimer linesearchTimer_;
  benchmark::RepeatedTimer computeControllerTimer_;
};

}  // namespace ocs2
