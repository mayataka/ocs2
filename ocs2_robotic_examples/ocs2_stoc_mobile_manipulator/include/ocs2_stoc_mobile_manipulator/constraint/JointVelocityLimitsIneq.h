#pragma once

#include <memory>

#include <ocs2_core/constraint/StateInputConstraint.h>

namespace ocs2 {
namespace mobile_manipulator {

class JointVelocityLimitsIneq final : public StateInputConstraint {
 public:
  explicit JointVelocityLimitsIneq(size_t inputDim, const vector_t& lowerBound, const vector_t& upperBound) 
    : StateInputConstraint(ConstraintOrder::Linear), inputDim_(inputDim), lowerBound_(lowerBound), upperBound_(upperBound) {}
  ~JointVelocityLimitsIneq() override = default;
  JointVelocityLimitsIneq* clone() const override { return new JointVelocityLimitsIneq(*this); }

  size_t getNumConstraints(scalar_t time) const override { return 2*inputDim_; }
  vector_t getValue(scalar_t time, const vector_t& state, const vector_t& input, const PreComputation&) const override;
  VectorFunctionLinearApproximation getLinearApproximation(scalar_t time, const vector_t& state, const vector_t& input,
                                                           const PreComputation&) const override;

 private:
  JointVelocityLimitsIneq(const JointVelocityLimitsIneq& other) = default;
  const size_t inputDim_;
  const vector_t lowerBound_, upperBound_;
};

}  // namespace mobile_manipulator
}  // namespace ocs2
