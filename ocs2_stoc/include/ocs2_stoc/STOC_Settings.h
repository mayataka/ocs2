#pragma once

#include <string>

#include <ocs2_core/Types.h>
#include <ocs2_stoc/riccati_recursion/RiccatiSolverMode.h>

namespace ocs2 {
namespace stoc {

struct Settings {
  // STOC settings
  size_t numIteration    = 10;       // Maximum number of Newton-type iterations
  scalar_t primalFeasTol = 1.0e-06;  // Termination condition : Primal feasibility, i.e., constraint violations, below this value
  scalar_t dualFeasTol   = 1.0e-06;  // Termination condition : Dual feasibility, i.e., KKT conditions, below this value
  scalar_t minPrimalStepSize = 0.0;  // Termination condition (failure) : Primal step size below thid value
  scalar_t minDualStepSize   = 0.0;  // Termination condition (failure) : Dual step size below thid value

  // Barrier strategy of the primal-dual interior point method. Conventions follows Ipopt.
  scalar_t initialBarrierParameter = 1.0e-03;       // Initial value of the barrier parameter
  scalar_t targetBarrierParameter  = 1.0e-03;       // Targer value of the barrier parameter. The barreir will decrease until reaches this value.
  scalar_t barrierLinearDecreaseFactor     = 0.2;   // Linear decrease factor of the barrier parameter, i.e., mu <- mu * factor.
  scalar_t barrierSuperlinearDecreasePower = 1.5;   // Superlinear decrease factor of the barrier parameter, i.e., mu <- mu ^ factor

  // If the current iterate satisfies the these two criteria, the barrier parameter is reduced (kkt-error based strategy of Ipopt).
  scalar_t barrierReductionPrimalFeasTol = 1.0e-03;  // Barrier reduction condition : Primal feasibility, i.e., constraint violations, below this value
  scalar_t barrierReductionDualFeasTol   = 1.0e-03;  // Barrier reduction condition : Dual feasibility, i.e., KKT conditions, below this value

  // Linesearch - step size rules
  scalar_t fractionToBoundaryMargin = 0.995;  // Margin of the fraction-to-boundary-rule for the step size selection 

  bool projectStateInputEqualityConstraints = true;  // Use a projection method to resolve the state-input constraint Cx+Du+e

  // controller type
  bool useFeedbackPolicy = true;  // true to use feedback, false to use feedforward

  // Discretization method
  scalar_t dt = 0.01;  // user-defined time discretization

  // STO strategy
  std::vector<std::pair<size_t, size_t>> stoEnabledModeSwitches; // Collection of the mode switches of which switching times will be optimized.
  scalar_t maxTimeInterval = 0.02; // Maximum time interval of the discretization in STO. 
  bool useOptimizedModeShceduleInReferenceManager = true; // If true, the optimized mode schedule is set to ReferenceManager after the optimization.
  scalar_t stoRegularization = 0.0;
  size_t skippedInitialStoModeSwitches = 0;
  size_t skippedFinalStoModeSwitches = 0;

  // If the current iterate satisfies the these two criteria, the mesh-refinement is performed.
  scalar_t meshRefinementPrimalFeasTol = 1.0e-02;  // Mesh refinement condition : Primal feasibility, i.e., constraint violations, below this value
  scalar_t meshRefinementDualFeasTol   = 1.0e-02;  // Mesh refinement condition : Dual feasibility, i.e., KKT conditions, below this value

  // Riccati options
  RiccatiSolverMode riccatiSolverMode = RiccatiSolverMode::Robust; // Robust: LDLT, Speed: LLT
  bool enableSwitchingTimeTrustRegion = true;                      // Use trust region-like regularization in the Riccati recursion for the switching times
  scalar_t switchingTimeTrustRegionRadius = 0.1;                   // Radius of the trust-region-like regularizaton of the switching times

  // Printing
  bool printSolverStatus = false;               // Print solver status 
  bool printSolverStatistics = false;           // Print benchmarking  
  bool printLinesearch = false;                 // Print linesearch (step-size) information
  bool printSwitchingTimeOptimization = false;  // Print switching time optimization-related information
  bool printDebugInfo = false;

  // Threading
  size_t nThreads = 4;
  int threadPriority = 99;
};

/**
 * Loads the STOC settings from a given file.
 *
 * @param [in] filename: File name which contains the configuration data.
 * @param [in] fieldName: Field name which contains the configuration data.
 * @param [in] verbose: Flag to determine whether to print out the loaded settings or not.
 * @return The settings
 */
Settings loadSettings(const std::string& filename, const std::string& fieldName = "stoc", bool verbose = true);

}  // namespace stoc
}  // namespace ocs2

