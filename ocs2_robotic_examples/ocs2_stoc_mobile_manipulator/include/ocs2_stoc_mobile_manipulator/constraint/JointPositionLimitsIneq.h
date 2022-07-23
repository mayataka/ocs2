#pragma once

#include <memory>

#include <ocs2_core/constraint/StateConstraint.h>

namespace ocs2 {
namespace mobile_manipulator {

class JointPositionLimitsIneq final : public StateConstraint {
 public:
  explicit JointPositionLimitsIneq(size_t armDim, const vector_t& lowerBound, const vector_t& upperBound) 
    : StateConstraint(ConstraintOrder::Linear), armDim_(armDim), lowerBound_(lowerBound), upperBound_(upperBound) {}
  ~JointPositionLimitsIneq() override = default;
  JointPositionLimitsIneq* clone() const override { return new JointPositionLimitsIneq(*this); }

  size_t getNumConstraints(scalar_t time) const override { return 2*armDim_; }
  vector_t getValue(scalar_t time, const vector_t& state, const PreComputation&) const override;
  VectorFunctionLinearApproximation getLinearApproximation(scalar_t time, const vector_t& state, const PreComputation&) const override;

 private:
  JointPositionLimitsIneq(const JointPositionLimitsIneq& other) = default;
  const size_t armDim_;
  const vector_t lowerBound_, upperBound_;
};

}  // namespace mobile_manipulator
}  // namespace ocs2
