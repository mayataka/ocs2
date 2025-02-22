#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_core/control/LinearController.h>
#include <ocs2_core/initialization/Initializer.h>
#include <ocs2_core/misc/Benchmark.h>
#include <ocs2_core/misc/LinearInterpolation.h>
#include <ocs2_core/misc/Numerics.h>
#include <ocs2_core/thread_support/ThreadPool.h>

#include <ocs2_oc/oc_solver/SolverBase.h>
#include <ocs2_oc/oc_data/DualSolution.h>
#include <ocs2_oc/oc_data/PerformanceIndex.h>
#include <ocs2_oc/oc_data/PrimalSolution.h>
#include <ocs2_oc/oc_data/ProblemMetrics.h>
#include <ocs2_oc/synchronized_module/ReferenceManagerInterface.h>
#include <ocs2_oc/synchronized_module/ReferenceManager.h>

#include <ocs2_ipm/model_data/ModelData.h>

#include <ocs2_ipm_oc/oc_problem/OptimalControlProblem.h>
#include <ocs2_ipm_oc/oc_data/PerformanceIndex.h>
#include <ocs2_ipm_oc/approximate_model/LinearQuadraticApproximator.h>
#include <ocs2_ipm_oc/approximate_model/LinearQuadraticDiscretizer.h>
#include <ocs2_ipm_oc/approximate_model/IpmVariablesEliminator.h>
#include <ocs2_ipm_oc/approximate_model/IpmVariablesRetriver.h>
#include <ocs2_ipm_oc/approximate_model/ConstraintProjection.h>

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

  /**
   * Sets the internal ReferenceManager for the switching time optimization. The switching time optimization is enabled only if the internal
   * ReferenceManager is set via this method.
   */
  void setInternalReferenceManager(std::shared_ptr<ReferenceManagerInterface> internalReferenceManagerPtr) {
    if (internalReferenceManagerPtr.get() == &getReferenceManager()) {
      throw std::runtime_error("[STOC] Resource of internalReferenceManagerPtr should be different from this->getReferenceManager()!");
    }
    internalReferenceManagerPtr_ = std::move(internalReferenceManagerPtr);
  }

  void reset() override;

  const OptimalControlProblem& getOptimalControlProblem() const override { return ocs2OcpDefinition_; }

  const PerformanceIndex& getPerformanceIndeces() const override { return getIterationsLog().back(); };

  size_t getNumIterations() const override { return totalNumIterations_; }

  const std::vector<PerformanceIndex>& getIterationsLog() const override;

  scalar_t getFinalTime() const override { return primalData_.primalSolution.timeTrajectory_.back(); };
  
  void getPrimalSolution(scalar_t finalTime, PrimalSolution* primalSolutionPtr) const override { *primalSolutionPtr = primalData_.primalSolution; }

  const ocs2::DualSolution& getDualSolution() const override {
    throw std::runtime_error("[STOC] getDualSolution() not available yet.");
  }

  const ProblemMetrics& getSolutionMetrics() const override {
    throw std::runtime_error("[STOC] getSolutionMetrics() not available yet.");
  }

  ScalarFunctionQuadraticApproximation getValueFunction(scalar_t time, const vector_t& state) const override {
    throw std::runtime_error("[STOC] getValueFunction() not available yet.");
  }

  ScalarFunctionQuadraticApproximation getHamiltonian(scalar_t time, const vector_t& state, const vector_t& input) override {
    throw std::runtime_error("[STOC] getHamiltonian() not available yet.");
  }

  vector_t getStateInputEqualityConstraintLagrangian(scalar_t time, const vector_t& state) const override {
    throw std::runtime_error("[STOC] getStateInputEqualityConstraintLagrangian() not available yet.");
  }

  MultiplierCollection getIntermediateDualSolution(scalar_t time) const override {
    throw std::runtime_error("[STOC] getIntermediateDualSolution() not available yet.");
  }

  const ipm::PerformanceIndex& getIpmPerformanceIndeces() const { return getIpmIterationsLog().back(); };

  const std::vector<ipm::PerformanceIndex>& getIpmIterationsLog() const;

  /**
   * Gets the benchmarking information as string
   */
  std::string getBenchmarkingInformation() const;

  /**
   * Const access to stoc settings
   */
  const stoc::Settings& getSettings() const { return settings_; }

 private:
  void runImpl(scalar_t initTime, const vector_t& initState, scalar_t finalTime) override;

  void runImpl(scalar_t initTime, const vector_t& initState, scalar_t finalTime, const ControllerBase* externalControllerPtr) override {
    if (externalControllerPtr == nullptr) {
      runImpl(initTime, initState, finalTime);
    } else {
      throw std::runtime_error("[STOC::run] This solver does not support external controller!");
    }
  }

  void runImpl(scalar_t initTime, const vector_t& initState, scalar_t finalTime, const PrimalSolution& primalSolution) override {
    // Copy all except the controller
    primalData_.primalSolution.timeTrajectory_ = primalSolution.timeTrajectory_;
    primalData_.primalSolution.stateTrajectory_ = primalSolution.stateTrajectory_;
    primalData_.primalSolution.inputTrajectory_ = primalSolution.inputTrajectory_;
    primalData_.primalSolution.postEventIndices_ = primalSolution.postEventIndices_;
    primalData_.primalSolution.modeSchedule_ = primalSolution.modeSchedule_;
    runImpl(initTime, initState, finalTime);
  }

  /** Run a task in parallel with settings.nThreads */
  void runParallel(std::function<void(int)> taskFunction);

  /** Initializes for the state-input trajectories */
  void initializeStateInputTrajectories(const std::vector<Grid>& timeDiscretization, const vector_t& initState,  
                                        vector_array_t& stateTrajectory, vector_array_t& inputTrajectory);

  /** Initializes for the costate trajectories */
  void initializeCostateTrajectories(const std::vector<Grid>& timeDiscretization, const vector_array_t& stateTrajectory, 
                                     const vector_array_t& inputTrajectory, vector_array_t& costateTrajectory);

  /** Creates quadratic approximation of the optimal control problem (OCP) around t, x(t), u(t). */
  ipm::PerformanceIndex approximateOptimalControlProblem(const std::vector<Grid>& timeDiscretization, const vector_t& initState, 
                                                         const vector_array_t& stateTrajectory, const vector_array_t& inputTrajectory,
                                                         const vector_array_t& costateTrajectory, 
                                                         std::vector<ipm::IpmVariables>& ipmVariablesTrajectory,
                                                         scalar_t barrierParameter, bool initIpmVariablesTrajectory);

  /** Creates quadratic approximation of the switching time optimization (STO) problem around current mode schedule. */
  ipm::PerformanceIndex approximateSwitchingTimeOptimizationProblem(scalar_t initTime, scalar_t finalTime, 
                                                                    const ModeSchedule& referenceModeSchedule, 
                                                                    const ModeSchedule& stoModeSchedule, 
                                                                    ipm::IpmVariables& stoIpmVariables,
                                                                    scalar_t barrierParameter, bool initStoIpmVariables);

  /** Summarizes the OCP approximation and STO approximatioin. */
  void summarizeModelData(const std::vector<Grid>& timeDiscretization);

  struct StepSizes {
    scalar_t primalStepSize, dualStepSize;
  };

  /** Computes the IPM directions and returns the primal and dual step sizes from OCP-data */
  StepSizes selectStepSizes(const std::vector<Grid>& timeDiscretization, const std::vector<ipm::IpmVariables>& ipmVariablesTrajectory,
                            const vector_array_t& dx, const vector_array_t& du, const vector_array_t& dlmd,
                            std::vector<ipm::IpmVariablesDirection>& ipmVariablesDirectionTrajectory, scalar_t fractionToBoundaryMargin); 

  /** Computes the IPM directions and returns the primal and dual step sizes from STO-data */
  StepSizes selectStepSizes(const ipm::IpmVariables& stoIpmVariables, const scalar_array_t& dts, 
                            ipm::IpmVariablesDirection& stoIpmVariablesDirection, scalar_t fractionToBoundaryMargin); 

  /** Updates the OCP-related variables */
  static void updateIterate(vector_array_t& x, vector_array_t& u, vector_array_t& lmd, std::vector<ipm::IpmVariables>& ipmVariablesTrajectory,
                            const vector_array_t& dx, const vector_array_t& du, const vector_array_t& dlmd, 
                            const std::vector<ipm::IpmVariablesDirection>& ipmVariablesDirectionTrajectory,
                            scalar_t primalStepSize, scalar_t dualStepSize);

  /** Updates the STO-related variables */
  void updateIterate(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, ModeSchedule& modeSchedule, 
                     ipm::IpmVariables& stoIpmVariables, const scalar_array_t& dts, 
                     const ipm::IpmVariablesDirection& stoIpmVariablesDirection, scalar_t primalStepSize, scalar_t dualStepSize);

  /** Set up the primal solution based on the optimized state and input trajectories */
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

  /** Determine convergence after a step */
  Convergence checkConvergence(size_t iteration, scalar_t barrierParameter, scalar_t maxTimeInterval, 
                               const ipm::PerformanceIndex& performanceIndex, scalar_t primalStepSize, scalar_t dualStepSize) const;

  /** Updates the barrier parameter */
  scalar_t updateBarrierParameter(scalar_t currentBarrierParameter,
                                  const ipm::PerformanceIndex& performanceIndex) const;

  // Problem definition
  stoc::Settings settings_;
  std::vector<ipm::OptimalControlProblem> optimalControlProblemStock_;
  OptimalControlProblem ocs2OcpDefinition_;
  std::unique_ptr<Initializer> initializerPtr_;

  // Data
  stoc::PrimalDataContainer primalData_;
  stoc::IpmDataContainer ipmData_;
  stoc::StoDataContainer stoData_;

  // riccati Recurision
  stoc::RiccatiRecursion riccatiRecursion_;

  // internal reference manager to manage mode schedules for the switching time optimization
  std::shared_ptr<ReferenceManagerInterface> internalReferenceManagerPtr_;

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
