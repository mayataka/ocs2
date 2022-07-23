#include <ocs2_stoc_mobile_manipulator/constraint/JointVelocityLimitsIneq.h>

namespace ocs2 {
namespace mobile_manipulator {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t JointVelocityLimitsIneq::getValue(scalar_t time, const vector_t& state, const vector_t& input, const PreComputation&) const {
  return input;
  vector_t constraint(2*inputDim_);
  constraint.head(inputDim_) = lowerBound_ - input;
  constraint.tail(inputDim_) = input - upperBound_; 
  return constraint;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximation JointVelocityLimitsIneq::getLinearApproximation(scalar_t time, const vector_t& state, 
                                                                                  const vector_t& input, 
                                                                                  const PreComputation& preComp) const {
  VectorFunctionLinearApproximation limits(2*inputDim_, state.rows(), inputDim_);
  limits.f = getValue(time, state, input, preComp);
  limits.dfdx.setZero();
  limits.dfdu.topRows(inputDim_) = - matrix_t::Identity(inputDim_, inputDim_);
  limits.dfdu.bottomRows(inputDim_) = matrix_t::Identity(inputDim_, inputDim_);
  return limits;
}

}  // namespace mobile_manipulator
}  // namespace ocs2
