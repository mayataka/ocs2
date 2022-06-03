#include <ocs2_ipm/approximate_model/LinearQuadraticDiscretizer.h>
#include <ocs2_ipm/ipm/InteriorPointMethod.h>

namespace ocs2 {
namespace ipm {

void discretizeIntermediateLQ(const scalar_t dt, const vector_t& state, const vector_t& state_next,
                              const DualVariable& dual, const DualVariable& dual_next, ModelData& modelData) {
  // hamiltonian
  modelData.hamiltonian.h = modelData.cost.f;
  modelData.hamiltonian.dhdx = modelData.cost.dfdx;
  modelData.hamiltonian.dhdu = modelData.cost.dfdu;
  modelData.hamiltonian.dhdx.noalias() += modelData.dynamics.dfdx.transpose() * dual_next.costate;
  modelData.hamiltonian.dhdu.noalias() += modelData.dynamics.dfdu.transpose() * dual_next.costate;
  modelData.hamiltonian.dfdt = modelData.dynamics.f;

  // dynamics 
  modelData.dynamics *= dt; 
  modelData.dynamics.f.noalias() += state;
  modelData.dynamics.f.noalias() -= state_next;

  // cost
  modelData.cost *= dt;
  modelData.cost.dfdx.noalias() += dual_next.costate;
  modelData.cost.dfdx.noalias() -= dual.costate;

  // state inequality constraints
  evalPerturbedResidual(modelData.stateIneqConstraint, dual.slackDualStateIneqConstraint, 
                        modelData.ipmDataStateIneqConstraint, modelData.cost, true, false);
  condenseSlackDual(modelData.stateIneqConstraint, dual.slackDualStateIneqConstraint, 
                    modelData.ipmDataStateIneqConstraint, modelData.cost, true, false);

  // state-input inequality constraints
  evalPerturbedResidual(modelData.stateInputIneqConstraint, dual.slackDualStateInputIneqConstraint, 
                        modelData.ipmDataStateInputIneqConstraint, modelData.cost, true, true);
  condenseSlackDual(modelData.stateInputIneqConstraint, dual.slackDualStateInputIneqConstraint, 
                    modelData.ipmDataStateInputIneqConstraint, modelData.cost, true, true);
}


void discretizePreJumpLQ(const vector_t& state, const vector_t& state_next,
                         const DualVariable& dual, const DualVariable& dual_next, ModelData& modelData) {
  // hamiltonian
  modelData.hamiltonian.h = modelData.cost.f;
  modelData.hamiltonian.dhdx = modelData.cost.dfdx;
  modelData.hamiltonian.dhdx.noalias() += modelData.dynamics.dfdx * dual_next.costate;
  modelData.hamiltonian.dfdt = modelData.dynamics.f;

  // dynamics 
  modelData.dynamics.f.noalias() += state;
  modelData.dynamics.f.noalias() -= state_next;

  // cost
  modelData.cost.dfdx.noalias() += dual_next.costate;
  modelData.cost.dfdx.noalias() -= dual.costate;

  // state inequality constraints
  evalPerturbedResidual(modelData.stateIneqConstraint, dual.slackDualStateIneqConstraint, 
                        modelData.ipmDataStateIneqConstraint, modelData.cost, true, false);
  condenseSlackDual(modelData.stateIneqConstraint, dual.slackDualStateIneqConstraint, 
                    modelData.ipmDataStateIneqConstraint, modelData.cost, true, false);
}


void discretizeFinalLQ(const DualVariable& dual, ModelData& modelData) {
  // hamiltonian
  modelData.hamiltonian.h = modelData.cost.f;
  modelData.hamiltonian.dhdx = modelData.cost.dfdx;

  // cost
  modelData.cost.dfdx.noalias() -= dual.costate;

  // state inequality constraints
  evalPerturbedResidual(modelData.stateIneqConstraint, dual.slackDualStateIneqConstraint, 
                        modelData.ipmDataStateIneqConstraint, modelData.cost, true, false);
  condenseSlackDual(modelData.stateIneqConstraint, dual.slackDualStateIneqConstraint, 
                    modelData.ipmDataStateIneqConstraint, modelData.cost, true, false);
}


void computeIntermediateDiscretizedMetrics(const scalar_t dt, const vector_t& state, const vector_t& state_next, 
                                           const DualVariable& dual, const ModelData& modelData, Metrics& metrics) {
  metrics.stateEquation = state + dt * modelData.dynamics.f - state_next;
  metrics.slackBarrier = slackLogBarrier(dual.slackDualStateIneqConstraint) + slackLogBarrier(dual.slackDualStateInputIneqConstraint);
}


void computePreJumpDiscretizedMetrics(const vector_t& state, const vector_t& state_next, const DualVariable& dual, 
                                      const ModelData& modelData, Metrics& metrics) {
  metrics.stateEquation = state + modelData.dynamics.f - state_next;
  metrics.slackBarrier = slackLogBarrier(dual.slackDualStateIneqConstraint);
}


void computeFinalDiscretizedMetrics(const vector_t& state, const DualVariable& dual, const ModelData& modelData, Metrics& metrics) { 
  metrics.slackBarrier = slackLogBarrier(dual.slackDualStateIneqConstraint);
}

}  // namespace ipm
}  // namespace ocs2
