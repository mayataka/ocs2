#include <ocs2_ipm/ipm/InteriorPointMethodSolution.h>

namespace ocs2 {

void setBarrier(scalar_t barrier, InteriorPointMethodSolution& solution) {
  for (auto& e : solution.slackDualTrajectory_) {
    setBarrier(barrier, e);
  }
}


std::string checkPositive(const InteriorPointMethodSolution& solution, const std::string& solutionName) {
  std::stringstream errorDescription;

  for (int i=0; i<solution.slackDualTrajectory_.size(); ++i) {
    errorDescription << checkPositive(solution.slackDualTrajectory_[i], 
                                      solutionName+"["+std::to_string(i)+"]");
  }
  return errorDescription.str();
}

}  // namespace ocs2
