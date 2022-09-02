#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/Hamiltonian.h>

namespace ocs2 {
namespace ipm {

/** Applies the change of input variables to a Hamiltonian */
void changeOfInputVariables(Hamiltonian& hamiltonian, const matrix_t& Pu, const matrix_t& Px = matrix_t(),
                            const vector_t& u0 = vector_t());

}  // namespace ipm
}  // namespace ocs2