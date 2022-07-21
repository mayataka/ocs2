#include "ocs2_stoc/STOC_Settings.h"

#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <ocs2_core/misc/LoadData.h>

namespace ocs2 {
namespace stoc {

Settings loadSettings(const std::string& filename, const std::string& fieldName, bool verbose) {
  boost::property_tree::ptree pt;
  boost::property_tree::read_info(filename, pt);

  Settings settings;

  if (verbose) {
    std::cerr << "\n #### STOC Settings:";
    std::cerr << "\n #### =============================================================================\n";
  }

  loadData::loadPtreeValue(pt, settings.numIteration, fieldName + ".numIteration", verbose);
  loadData::loadPtreeValue(pt, settings.primalFeasTol, fieldName + ".primalFeasTol", verbose);
  loadData::loadPtreeValue(pt, settings.dualFeasTol, fieldName + ".dualFeasTol", verbose);
  loadData::loadPtreeValue(pt, settings.initialBarrierParameter, fieldName + ".initialBarrierParameter", verbose);
  loadData::loadPtreeValue(pt, settings.targetBarrierParameter, fieldName + ".targetBarrierParameter", verbose);
  loadData::loadPtreeValue(pt, settings.barrierLinearDecreaseFactor, fieldName + ".barrierLinearDecreaseFactor", verbose);
  loadData::loadPtreeValue(pt, settings.barrierSuperlinearDecreasePower, fieldName + ".barrierSuperlinearDecreasePower", verbose);
  loadData::loadPtreeValue(pt, settings.barrierReductionPrimalFeasTol, fieldName + ".barrierReductionPrimalFeasTol", verbose);
  loadData::loadPtreeValue(pt, settings.barrierReductionDualFeasTol, fieldName + ".barrierReductionDualFeasTol", verbose);
  loadData::loadPtreeValue(pt, settings.fractionToBoundaryMargin, fieldName + ".fractionToBoundaryMargin", verbose);
  loadData::loadPtreeValue(pt, settings.useFeedbackPolicy, fieldName + ".useFeedbackPolicy", verbose);
  loadData::loadPtreeValue(pt, settings.dt, fieldName + ".dt", verbose);
  loadData::loadPtreeValue(pt, settings.printSolverStatus, fieldName + ".printSolverStatus", verbose);
  loadData::loadPtreeValue(pt, settings.printSolverStatistics, fieldName + ".printSolverStatistics", verbose);
  loadData::loadPtreeValue(pt, settings.printLinesearch, fieldName + ".printLinesearch", verbose);
  loadData::loadPtreeValue(pt, settings.nThreads, fieldName + ".nThreads", verbose);
  loadData::loadPtreeValue(pt, settings.threadPriority, fieldName + ".threadPriority", verbose);

  if (verbose) {
    std::cerr << " #### =============================================================================" << std::endl;
  }

  return settings;
}
}  // namespace stoc
}  // namespace ocs2