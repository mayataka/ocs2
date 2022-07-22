#include "ocs2_ipm/approximate_model/ChangeOfInputVariables.h"

namespace ocs2 {
namespace ipm {

void changeOfInputVariables(Hamiltonian& hamiltonian, const matrix_t& Pu, const matrix_t& Px, const vector_t& u0) {
  const bool hasPx(Px.size() > 0);
  const bool hasu0(u0.size() > 0);

  if (hasPx) {
    hamiltonian.dhdx.noalias() += Px.transpose() * hamiltonian.dhdu; // Before adapting dhdu!
  }

  if (hasu0) {
    hamiltonian.h += u0.dot(hamiltonian.dhdu); // Before adapting dfdu!
  }

  hamiltonian.dhdu = Pu.transpose() * hamiltonian.dhdu;  
}

}  // namespace ipm
}  // namespace ocs2
