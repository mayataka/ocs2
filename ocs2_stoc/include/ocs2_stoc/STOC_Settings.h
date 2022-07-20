#pragma once

#include <string>

#include <ocs2_core/Types.h>

namespace ocs2 {
namespace stoc {

struct Settings {
  // STOC settings
  size_t numIteration    = 10;  // Maximum number of Newton-type iterations
  scalar_t primalFeasTol = 1e-6; // Termination condition : Primal feasibility, i.e., constraint violations, below this value
  scalar_t dualFeasTol   = 1e-4; // Termination condition : Dual feasibility, i.e., KKT conditions, below this value
  scalar_t minPrimalStepSize = 1.0e-03;
  scalar_t minDualStepSize = 1.0e-03;

  // Barrier parameter for primal-dual interior point method
  scalar_t initialBarrierParameter = 1.0e-01;
  scalar_t targetBarrierParameter = 1.0e-03;
  scalar_t barrierReductionRate = 0.5;

  // Linesearch - step size rules
  scalar_t fractionToBoundaryMargin = 0.995;  // Margin of the fraction-to-boundary-rule for the step size selection 

  // controller type
  bool useFeedbackPolicy = true;  // true to use feedback, false to use feedforward

  // Discretization method
  scalar_t dt = 0.01;  // user-defined time discretization

  // Printing
  bool printSolverStatus = false;      // Print HPIPM status after solving the QP subproblem
  bool printSolverStatistics = false;  // Print benchmarking of the multiple shooting method
  bool printLinesearch = false;        // Print linesearch information

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

