#include <ocs2_stoc_mobile_manipulator/constraint/JointPositionLimitsIneq.h>

namespace ocs2 {
namespace mobile_manipulator {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t JointPositionLimitsIneq::getValue(scalar_t time, const vector_t& state, const PreComputation&) const {
  vector_t constraint(2*armDim_);
  constraint.head(armDim_) = lowerBound_ - state.tail(armDim_);
  constraint.tail(armDim_) = state.tail(armDim_) - upperBound_; 
  return constraint;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximation JointPositionLimitsIneq::getLinearApproximation(scalar_t time, const vector_t& state,
                                                                                  const PreComputation& preComp) const {
  const size_t baseDim = state.rows() - armDim_;

  VectorFunctionLinearApproximation limits(2*armDim_, state.rows(), 0);
  limits.f = getValue(time, state, preComp);
  limits.dfdx.topLeftCorner(2*armDim_, baseDim).setZero();
  limits.dfdx.topRightCorner(2*armDim_, armDim_).topRows(armDim_) = - matrix_t::Identity(armDim_, armDim_);
  limits.dfdx.topRightCorner(2*armDim_, armDim_).bottomRows(armDim_) = matrix_t::Identity(armDim_, armDim_);
  return limits;
}

}  // namespace mobile_manipulator
}  // namespace ocs2
