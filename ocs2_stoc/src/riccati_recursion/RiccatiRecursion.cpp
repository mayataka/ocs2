#include <ocs2_stoc/riccati_recursion/RiccatiRecursion.h>

#include <cassert>

namespace ocs2 {
namespace stoc {

RiccatiRecursion::RiccatiRecursion(scalar_t dts0_max)
  : riccati_(),
    lqrPolicy_(),
    stoPolicy_(),
    backwardRecursion_(dts0_max) {
}


void RiccatiRecursion::resize(size_t N) {
  if (riccati_.size() < N+1) {
    riccati_.resize(N+1);
    lqrPolicy_.resize(N);
    stoPolicy_.resize(N);
  }
}


void RiccatiRecursion::backwardRecursion(const std::vector<Grid>& timeDiscretizationGrid, std::vector<ipm::ModelData>& modelData) {
  const size_t N = timeDiscretizationGrid.size() - 1;
  resize(N);
  backwardRecursion_.computeFinal(modelData[N], riccati_[N]);
  for (int i=N-1; i>=0; --i) {
    const bool sto     = timeDiscretizationGrid[i].sto;
    const bool stoNext = timeDiscretizationGrid[i].stoNext;
    if (timeDiscretizationGrid[i].event == Grid::Event::PreEvent) {
      assert(timeDiscretizationGrid[i+1].event == Grid::Event::PostEvent);
      backwardRecursion_.computeIntermediate(riccati_[i+1], modelData[i], riccati_[i], lqrPolicy_[i], sto, stoNext);
    }
    else {
      const auto phase = timeDiscretizationGrid[i].phase;
      backwardRecursion_.computePreJump(riccati_[i+1], modelData[i], riccati_[i], lqrPolicy_[i], stoPolicy_[phase+1], sto, stoNext);
    }
  }
}


void RiccatiRecursion::forwardRecursion(const std::vector<Grid>& timeDiscretizationGrid, const std::vector<ipm::ModelData>& modelData,
                                        vector_array_t& stateTrajectory, vector_array_t& inputTrajectory, 
                                        vector_array_t& costateTrajectory, scalar_array_t& switchingTimes) {
  const size_t N = timeDiscretizationGrid.size() - 1;
  const auto phase = timeDiscretizationGrid[0].phase;
  const bool sto   = timeDiscretizationGrid[0].sto;
  switchingTimes[phase] = 0.0;
  if (sto) {
    // This is the switching time from phase to phase+1
    switchingTimes[phase+1] = ForwardRiccatiRecursion::computeSwitchingTime(stoPolicy_[phase], stateTrajectory[0], 
                                                                            switchingTimes[phase], false);
  }
  for (int i=0; i<N; ++i) {
    const auto phase     = timeDiscretizationGrid[i].phase;
    const auto phaseNext = timeDiscretizationGrid[i+1].phase;
    const bool sto     = timeDiscretizationGrid[i].sto;
    const bool stoNext = timeDiscretizationGrid[i].stoNext;
    ForwardRiccatiRecursion::computeInput(lqrPolicy_[i], stateTrajectory[i], inputTrajectory[i], 
                                          switchingTimes[phase], switchingTimes[phase+1], sto, stoNext);
    ForwardRiccatiRecursion::computeState(modelData[i], stateTrajectory[i], inputTrajectory[i], stateTrajectory[i+1], 
                                          switchingTimes[phase], switchingTimes[phase+1], sto);
    ForwardRiccatiRecursion::computeCostate(riccati_[i], stateTrajectory[i], costateTrajectory[i],
                                            switchingTimes[phase], switchingTimes[phase+1], sto, stoNext);
    if (phase != phaseNext) {
      assert(phase+1 == phaseNext);
      // In this case, the new mode becomes active at the grid i+1. 
      // This switch is indexed by the phase.
      switchingTimes[phase+1] = ForwardRiccatiRecursion::computeSwitchingTime(stoPolicy_[phase+1], stateTrajectory[i+1], 
                                                                              switchingTimes[phase], stoNext);
    }
  }
  ForwardRiccatiRecursion::computeCostate(riccati_[N], stateTrajectory[N], costateTrajectory[N]);
}

} // namespace stoc
} // namespace ocs2