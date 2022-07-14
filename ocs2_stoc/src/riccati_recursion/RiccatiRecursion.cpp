#include <ocs2_stoc/riccati_recursion/RiccatiRecursion.h>

#include <cassert>

namespace ocs2 {
namespace stoc {

RiccatiRecursion::RiccatiRecursion(const size_t nx, const size_t nu, const size_t N,
                                   const size_t maxNumEvents, const scalar_t dts0_max)
  : N_(N),
    riccati_(),
    riccatiPreEvent_(),
    lqrPolicy_(),
    stoPolicy_(),
    backwardRecursion_(nx, nu, dts0_max) {
  if (nx <= 0) {
    throw std::out_of_range("nx must be positive!");
  }
  if (nu < 0) {
    throw std::out_of_range("nu must be non-negative!");
  }
  if (N <= 0) {
    throw std::out_of_range("N must be positive!");
  }
  if (maxNumEvents < 0) {
    throw std::out_of_range("maxNumEvents must be non-negative!");
  }
  if (dts0_max < 0) {
    throw std::out_of_range("dts0_max must be non-negative!");
  }
  riccati_.resize(N+1);
  for (auto& e : riccati_) {
    e.resize(nx, nu);
    e.setZero();
  }
  riccatiPreEvent_.resize(maxNumEvents+1);
  for (int i=0; i<maxNumEvents+1; ++i) {
    riccatiPreEvent_[i].resize(nx, nu);
    riccatiPreEvent_[i].setZero();
  }
  lqrPolicy_.resize(N);
  for (auto& e : lqrPolicy_) {
    e.resize(nx, nu);
    e.setZero();
  }
  stoPolicy_.resize(maxNumEvents+1);
  for (int i=0; i<maxNumEvents; ++i) {
    stoPolicy_[i].resize(nx);
    stoPolicy_[i].setZero();
  }
}


RiccatiRecursion::RiccatiRecursion() {
}


void RiccatiRecursion::backwardRecursion(const std::vector<Grid>& timeDiscretizationGrid, std::vector<ipm::ModelData>& modelData) {
  backwardRecursion_.compute(modelData[N_].cost, riccati_[N_]);
  for (int i=N_-1; i>=0; --i) {
    const auto phase     = timeDiscretizationGrid[i].phase;
    const auto phaseNext = timeDiscretizationGrid[i+1].phase;
    const bool sto       = timeDiscretizationGrid[i].sto;
    const bool stoNext   = timeDiscretizationGrid[i].stoNext;
    if (phase != phaseNext) {
      assert(phase+1 == phaseNext);
      // In this case, the new mode becomes active at the grid i+1. 
      // This switch is indexed by the phase.
      const bool stoNextNext = timeDiscretizationGrid[i].stoNextNext;
      backwardRecursion_.phaseTransition(riccati_[i+1], riccatiPreEvent_[phase+1], stoPolicy_[phase+1], stoNextNext);
      backwardRecursion_.compute(riccatiPreEvent_[phase+1], modelData[i].dynamics, modelData[i].cost, modelData[i].hamiltonian, 
                                 riccati_[i], lqrPolicy_[i], sto, stoNext);
    }
    else {
      backwardRecursion_.compute(riccati_[i+1], modelData[i].dynamics, modelData[i].cost, modelData[i].hamiltonian, 
                                 riccati_[i], lqrPolicy_[i], sto, stoNext);
    }
  }
  const auto phase = timeDiscretizationGrid[0].phase;
  const bool sto   = timeDiscretizationGrid[0].sto;
  if (sto) {
    const bool stoNext = timeDiscretizationGrid[0].stoNext;
    backwardRecursion_.phaseTransition(riccati_[0], riccatiPreEvent_[phase], stoPolicy_[phase], stoNext);
  }
}


void RiccatiRecursion::forwardRecursion(const std::vector<Grid>& timeDiscretizationGrid, const std::vector<ipm::ModelData>& modelData,
                                        vector_array_t& stateTrajectory, vector_array_t& inputTrajectory, 
                                        vector_array_t& costateTrajectory, scalar_array_t& switchingTimes) {
  const auto phase = timeDiscretizationGrid[0].phase;
  const bool sto   = timeDiscretizationGrid[0].sto;
  switchingTimes[phase] = 0.0;
  if (sto) {
    // This is the switching time from phase to phase+1
    switchingTimes[phase+1] = ForwardRiccatiRecursion::computeSwitchingTime(stoPolicy_[phase], stateTrajectory[0], 
                                                                            switchingTimes[phase], false);
  }
  for (int i=0; i<N_; ++i) {
    const auto phase     = timeDiscretizationGrid[i].phase;
    const auto phaseNext = timeDiscretizationGrid[i+1].phase;
    const bool sto     = timeDiscretizationGrid[i].sto;
    const bool stoNext = timeDiscretizationGrid[i].stoNext;
    ForwardRiccatiRecursion::computeInput(lqrPolicy_[i], stateTrajectory[i], inputTrajectory[i], 
                                          switchingTimes[phase], switchingTimes[phase+1], sto, stoNext);
    ForwardRiccatiRecursion::computeState(modelData[i].dynamics, modelData[i].hamiltonian, stateTrajectory[i], inputTrajectory[i], 
                                          stateTrajectory[i+1], switchingTimes[phase], switchingTimes[phase+1], sto);
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
  ForwardRiccatiRecursion::computeCostate(riccati_[N_], stateTrajectory[N_], costateTrajectory[N_]);
}

} // namespace stoc
} // namespace ocs2