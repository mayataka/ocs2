#include "ocs2_stoc/riccati_recursion/RiccatiRecursion.h"

#include <cassert>

namespace ocs2 {
namespace stoc {

RiccatiRecursion::RiccatiRecursion(RiccatiSolverMode riccatiSolverMode, scalar_t switchingTimeTrustRegion, 
                                   bool enableSwitchingTimeTrustRegion)
  : riccatiData_(),
    lqrPolicy_(),
    stoPolicy_(),
    backwardRecursion_(riccatiSolverMode, switchingTimeTrustRegion, enableSwitchingTimeTrustRegion) {
}


void RiccatiRecursion::resize(size_t N) {
  if (riccatiData_.size() < N+1) {
    riccatiData_.resize(N+1);
    lqrPolicy_.resize(N);
    stoPolicy_.resize(N);
  }
}


void RiccatiRecursion::backwardRecursion(const std::vector<Grid>& timeDiscretizationGrid, std::vector<ipm::ModelData>& modelData) {
  const size_t N = timeDiscretizationGrid.size() - 1;
  resize(N);
  backwardRecursion_.computeFinal(modelData[N], riccatiData_[N]);
  for (int i=N-1; i>=0; --i) {
    const bool sto     = timeDiscretizationGrid[i].sto;
    const bool stoNext = timeDiscretizationGrid[i].stoNext;
    if (timeDiscretizationGrid[i].event == Grid::Event::PreEvent) {
      const auto phase = timeDiscretizationGrid[i].phase;
      const bool stoNextNext = timeDiscretizationGrid[i].stoNextNext;
      backwardRecursion_.computePreJump(riccatiData_[i+1], modelData[i], riccatiData_[i], lqrPolicy_[i], stoPolicy_[phase+1], sto, stoNext, 
                                        stoNextNext);
    } else {
      backwardRecursion_.computeIntermediate(riccatiData_[i+1], modelData[i], riccatiData_[i], lqrPolicy_[i], sto, stoNext);
    }
  }
  if (timeDiscretizationGrid[0].sto) {
    const auto phase = timeDiscretizationGrid[0].phase;
    backwardRecursion_.modifyPreJump(riccatiData_[0], stoPolicy_[phase], timeDiscretizationGrid[0].stoNext);
    riccatiData_[0].Phi.setZero();
  }
}


void RiccatiRecursion::forwardRecursion(const std::vector<Grid>& timeDiscretizationGrid, const std::vector<ipm::ModelData>& modelData,
                                        vector_array_t& stateTrajectory, vector_array_t& inputTrajectory, 
                                        vector_array_t& costateTrajectory, scalar_array_t& switchingTimesOutput) {
  const size_t N = timeDiscretizationGrid.size() - 1;
  const size_t numPhases = timeDiscretizationGrid.back().phase + 1;
  scalar_array_t switchingTimes(numPhases, 0.0);
  switchingTimes[timeDiscretizationGrid[0].phase] = 0.0;
  if (timeDiscretizationGrid[0].sto) {
    // This is the switching time from phase to phase+1
    const auto phase = timeDiscretizationGrid[0].phase;
    constexpr bool stoPrev = false;
    switchingTimes[phase+1] = ForwardRiccatiRecursion::computeSwitchingTime(stoPolicy_[phase], stateTrajectory[0], 
                                                                            switchingTimes[phase], stoPrev);
  }
  for (int i=0; i<N; ++i) {
    const auto phase     = timeDiscretizationGrid[i].phase;
    const bool sto     = timeDiscretizationGrid[i].sto;
    const bool stoNext = timeDiscretizationGrid[i].stoNext;
    if (timeDiscretizationGrid[i].event == Grid::Event::PostEvent && stoNext) {
      switchingTimes[phase+1] = ForwardRiccatiRecursion::computeSwitchingTime(stoPolicy_[phase], stateTrajectory[i], 
                                                                              switchingTimes[phase], sto);
    }
    ForwardRiccatiRecursion::computeInput(lqrPolicy_[i], stateTrajectory[i], inputTrajectory[i], 
                                          switchingTimes[phase], switchingTimes[phase+1], sto, stoNext);
    ForwardRiccatiRecursion::computeState(modelData[i], stateTrajectory[i], inputTrajectory[i], stateTrajectory[i+1], 
                                          switchingTimes[phase], switchingTimes[phase+1], sto);
    ForwardRiccatiRecursion::computeCostate(riccatiData_[i], stateTrajectory[i], costateTrajectory[i],
                                            switchingTimes[phase], switchingTimes[phase+1], sto, stoNext);
  }
  ForwardRiccatiRecursion::computeCostate(riccatiData_[N], stateTrajectory[N], costateTrajectory[N]);
  switchingTimesOutput.clear(); switchingTimesOutput.resize(numPhases-1);
  for (int i=0; i<switchingTimesOutput.size(); ++i) {
    switchingTimesOutput[i] = switchingTimes[i+1];
  }
}

} // namespace stoc
} // namespace ocs2