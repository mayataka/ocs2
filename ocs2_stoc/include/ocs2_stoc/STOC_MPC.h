#pragma once

#include <ocs2_mpc/MPC_BASE.h>

#include <ocs2_stoc/STOC.h>

namespace ocs2 {
class STOC_MPC final : public MPC_BASE {
 public:
  /**
   * Constructor
   *
   * @param mpcSettings : settings for the mpc wrapping of the solver. Do not use this for maxIterations and stepsize, use settings directly.
   * @param settings : settings for the STOC solver.
   * @param [in] optimalControlProblem: The optimal control problem formulation.
   * @param [in] initializer: This class initializes the state-input for the time steps that no controller is available.
   */
  STOC_MPC(mpc::Settings mpcSettings, stoc::Settings settings, const ipm::OptimalControlProblem& optimalControlProblem,
           const Initializer& initializer)
      : MPC_BASE(std::move(mpcSettings)) {
    solverPtr_.reset(new STOC(std::move(settings), optimalControlProblem, initializer));
  };

  ~STOC_MPC() override = default;

  STOC* getSolverPtr() override { return solverPtr_.get(); }
  const STOC* getSolverPtr() const override { return solverPtr_.get(); }

 protected:
  void calculateController(scalar_t initTime, const vector_t& initState, scalar_t finalTime) override {
    if (settings().coldStart_) {
      solverPtr_->reset();
    }
    solverPtr_->run(initTime, initState, finalTime);
  }

 private:
  std::unique_ptr<STOC> solverPtr_;
};
}  // namespace ocs2
