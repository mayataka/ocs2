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


void RiccatiRecursion::backwardRecursion(const DiscreteTimeModeSchedule& modeSchedule, std::vector<ipm::ModelData>& modelData) {
  backwardRecursion_.compute(modelData[N_].cost, riccati_[N_]);
  for (int i=N_-1; i>=0; --i) {
    const auto phase     = modeSchedule.phaseAtTimeStage(i);
    const auto phaseNext = modeSchedule.phaseAtTimeStage(i+1);
    const bool sto     = modeSchedule.isStoEnabledAtPhase(phase);
    const bool stoNext = modeSchedule.isStoEnabledAtPhase(phase+1); 
    if (phase != phaseNext) {
      assert(phase+1 == phaseNext);
      // In this case, the new mode becomes active at the grid i+1. 
      // This switch is indexed by the phase.
      const bool stoNextNext = modeSchedule.isStoEnabledAtPhase(phase+2); 
      backwardRecursion_.phaseTransition(riccati_[i+1], riccatiPreEvent_[phase+1], stoPolicy_[phase+1], stoNextNext);
      backwardRecursion_.compute(riccatiPreEvent_[phase+1], modelData[i].dynamics, modelData[i].cost, modelData[i].hamiltonian, 
                                 riccati_[i], lqrPolicy_[i], sto, stoNext);
    }
    else {
      backwardRecursion_.compute(riccati_[i+1], modelData[i].dynamics, modelData[i].cost, modelData[i].hamiltonian, 
                                 riccati_[i], lqrPolicy_[i], sto, stoNext);
    }
  }
  const auto phase = modeSchedule.phaseAtTimeStage(0);
  if (modeSchedule.isStoEnabledAtPhase(phase)) {
    const bool stoNext = modeSchedule.isStoEnabledAtPhase(phase+1);
    backwardRecursion_.phaseTransition(riccati_[0], riccatiPreEvent_[phase], stoPolicy_[phase], stoNext);
  }
}


void RiccatiRecursion::forwardRecursion(const DiscreteTimeModeSchedule& modeSchedule, const std::vector<ipm::ModelData>& modelData, 
                                        vector_array_t& stateTrajectory, vector_array_t& inputTrajectory, 
                                        vector_array_t& costateTrajectory, scalar_array_t& switchingTimes) {
  const auto phase = modeSchedule.phaseAtTimeStage(0);
  switchingTimes[phase] = 0.0;
  if (modeSchedule.isStoEnabledAtPhase(phase)) {
    // This is the switching time from phase to phase+1
    switchingTimes[phase+1] = ForwardRiccatiRecursion::computeSwitchingTime(stoPolicy_[phase], stateTrajectory[0], 
                                                                            switchingTimes[phase], false);
  }
  for (int i=0; i<N_; ++i) {
    const auto phase     = modeSchedule.phaseAtTimeStage(i);
    const auto phaseNext = modeSchedule.phaseAtTimeStage(i+1);
    const bool sto         = modeSchedule.isStoEnabledAtPhase(phase);
    const bool stoNext     = modeSchedule.isStoEnabledAtPhase(phase+1);
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