#include <ocs2_ipm/approximate_model/LinearQuadraticDiscretizer.h>

namespace ocs2 {
namespace ipm {

void discretizeIntermediateLQ(const scalar_t dt, const vector_t& state, const vector_t& stateNext, const vector_t& costate, 
                              const vector_t& costateNext, ModelData& modelData) {
  // cost 
  modelData.cost.dfdx.noalias() += modelData.dynamics.dfdx.transpose() * costateNext; 
  modelData.cost.dfdu.noalias() += modelData.dynamics.dfdu.transpose() * costateNext; 

  // hamiltonian
  modelData.hamiltonian.h = modelData.cost.f;
  modelData.hamiltonian.dhdx = modelData.cost.dfdx;
  modelData.hamiltonian.dhdu = modelData.cost.dfdu;
  modelData.hamiltonian.dfdt = modelData.dynamics.f;

  // dynamics 
  modelData.dynamics *= dt; 
  modelData.dynamics.f.noalias() += state;
  modelData.dynamics.f.noalias() -= stateNext;

  // cost
  modelData.cost *= dt;
  modelData.cost.dfdx.noalias() += costateNext;
  modelData.cost.dfdx.noalias() -= costate;
}


void discretizePreJumpLQ(const vector_t& state, const vector_t& stateNext, const vector_t& costate, 
                         const vector_t& costateNext, ModelData& modelData) {
  // cost 
  modelData.cost.dfdx.noalias() += modelData.dynamics.dfdx.transpose() * costateNext; 

  // hamiltonian
  modelData.hamiltonian.h = modelData.cost.f;
  modelData.hamiltonian.dhdx = modelData.cost.dfdx;
  modelData.hamiltonian.dfdt = modelData.dynamics.f;

  // dynamics 
  modelData.dynamics.f.noalias() += state;
  modelData.dynamics.f.noalias() -= stateNext;

  // cost
  modelData.cost.dfdx.noalias() += costateNext;
  modelData.cost.dfdx.noalias() -= costate;
}


void discretizeFinalLQ(const vector_t& costate, ModelData& modelData) {
  // hamiltonian
  modelData.hamiltonian.h = modelData.cost.f;
  modelData.hamiltonian.dhdx = modelData.cost.dfdx;

  // cost
  modelData.cost.dfdx.noalias() -= costate;
}

}  // namespace ipm
}  // namespace ocs2
