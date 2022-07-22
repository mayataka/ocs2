#pragma once

#include <ocs2_core/misc/LinearInterpolation.h>
#include <ocs2_ipm/model_data/ModelData.h>

/*
 * @file
 * The linear interpolation of cost inpute-state derivative, Pm, at index-alpha pair given modelDataTrajectory can
 * be computed as:
 *
 * LinearInterpolation::interpolate(indexAlpha, Pm, &modelDataTrajectory, model_data::cost_dfdux);
 */

/*
 * Declares an access function of name FIELD such as time, dynamics, dynamicsBias, ...
 * For example the signature of function for dynamics is:
 * const vector_t& dynamics(const std::vector<ocs2::ModelData>& vec, size_t n) {
 *   return vec[n].dynamics;
 * }
 */
#define CREATE_INTERPOLATION_ACCESS_FUNCTION(FIELD) \
  inline auto FIELD(const std::vector<ocs2::ipm::ModelData>& vec, size_t ind)->const decltype(vec[ind].FIELD)& { return vec[ind].FIELD; }

#define CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(FIELD, SUBFIELD)                                                            \
  inline auto FIELD##_##SUBFIELD(const std::vector<ocs2::ipm::ModelData>& vec, size_t ind)->const decltype(vec[ind].FIELD.SUBFIELD)& { \
    return vec[ind].FIELD.SUBFIELD;                                                                                               \
  }

namespace ocs2 {
namespace ipm {
namespace model_data {

/**
 * Access method for different subfields of the ModelData.
 */

// Time
CREATE_INTERPOLATION_ACCESS_FUNCTION(time)

// Dynamics
CREATE_INTERPOLATION_ACCESS_FUNCTION(dynamicsBias)
CREATE_INTERPOLATION_ACCESS_FUNCTION(dynamicsCovariance)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(dynamics, f)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(dynamics, dfdx)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(dynamics, dfdu)

// Cost
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(cost, f)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(cost, dfdx)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(cost, dfdxx)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(cost, dfdu)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(cost, dfduu)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(cost, dfdux)

// State equality constraints
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(stateEqConstraint, f)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(stateEqConstraint, dfdx)

// State-input equality constraints
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(stateInputEqConstraint, f)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(stateInputEqConstraint, dfdx)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(stateInputEqConstraint, dfdu)

// State inequality constraints
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(stateIneqConstraint, f)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(stateIneqConstraint, dfdx)

// State-input inequality constraints
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(stateInputIneqConstraint, f)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(stateInputIneqConstraint, dfdx)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(stateInputIneqConstraint, dfdu)

// Hamiltonian
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(hamiltonian, h)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(hamiltonian, dhdt)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(hamiltonian, dhdx)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(hamiltonian, dhdu)
CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD(hamiltonian, dfdt)

}  // namespace model_data
}  // namespace ipm
}  // namespace ocs2

#undef CREATE_INTERPOLATION_ACCESS_FUNCTION
#undef CREATE_INTERPOLATION_ACCESS_FUNCTION_SUBFIELD
