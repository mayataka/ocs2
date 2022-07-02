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
#include <ocs2_ipm/model_data/ModelData.h>
// #include <ocs2_ipm/model_data/ModelDataLinearInterpolation.h> :TODO
#include <ocs2_ipm/approximate_model/LinearQuadraticApproximator.h>
#include <ocs2_ipm/approximate_model/LinearQuadraticDiscretizer.h>
#include <ocs2_ipm/approximate_model/DualVariableReconstructor.h>

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

  scalar_t getFinalTime() const override;

  void getPrimalSolution(scalar_t finalTime, PrimalSolution* primalSolutionPtr) const override;

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

  PerformanceIndex approximateOptimalControlProblem(const vector_t& initState, const std::vector<AnnotatedTime>& timeDiscretization);

  std::pair<scalar_t, scalar_t> fractionToBoundaryRule(const std::vector<AnnotatedTime>& timeDiscretization, 
                                                       const vector_array_t& dx, const vector_array_t& du, const scalar_array_t& dts); 

  // Problem definition
  stoc::Settings settings_;
  std::vector<ipm::OptimalControlProblem> optimalControlProblemStock_;
  std::unique_ptr<Initializer> initializerPtr_;

  // Data
  stoc::PrimalDataContainer primalData_;
  stoc::DualDataContainer dualData_;

  // riccati Recurision
  stoc::RiccatiRecursion riccatiRecursion_;

  // Directions
  vector_array_t dx_;
  vector_array_t du_;
  vector_array_t dlmd_;
  scalar_array_t dts_;
  std::vector<ipm::DualVariableDirection> dualDirectionTrajectory_;

  // Threading
  ThreadPool threadPool_;

  // Iteration performance log
  std::vector<PerformanceIndex> performanceIndeces_;

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
