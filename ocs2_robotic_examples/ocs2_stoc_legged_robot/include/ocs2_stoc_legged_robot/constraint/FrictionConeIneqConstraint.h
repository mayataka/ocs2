#pragma once

#include <ocs2_centroidal_model/CentroidalModelInfo.h>
#include <ocs2_core/constraint/StateInputConstraint.h>

#include "ocs2_legged_robot/common/Types.h"
#include "ocs2_legged_robot/reference_manager/SwitchedModelReferenceManager.h"
#include "ocs2_legged_robot/constraint/FrictionConeConstraint.h"

namespace ocs2 {
namespace legged_robot {

/**
 * Implements the hard constraint h(t,x,u) >= 0
 *
 * frictionCoefficient * (Fz + gripperForce) - sqrt(Fx * Fx + Fy * Fy + regularization) >= 0
 *
 * The gripper force shifts the origin of the friction cone down in z-direction by the amount of gripping force available. This makes it
 * possible to produce tangential forces without applying a regular normal force on that foot, or to "pull" on the foot with magnitude up to
 * the gripping force.
 *
 * The regularization prevents the constraint gradient / hessian to go to infinity when Fx = Fz = 0. It also creates a parabolic safety
 * margin to the friction cone. For example: when Fx = Fy = 0 the constraint zero-crossing will be at Fz = 1/frictionCoefficient *
 * sqrt(regularization) instead of Fz = 0
 *
 */
class FrictionConeIneqConstraint final : public StateInputConstraint {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  using Config = FrictionConeConstraint::Config;

  /**
   * Constructor
   * @param [in] referenceManager : Switched model ReferenceManager.
   * @param [in] config : Friction model settings.
   * @param [in] contactPointIndex : The 3 DoF contact index.
   * @param [in] info : The centroidal model information.
   */
  FrictionConeIneqConstraint(const SwitchedModelReferenceManager& referenceManager, Config config, size_t contactPointIndex,
                             CentroidalModelInfo info);

  // FrictionConeIneqConstraint(const FrictionConeConstraint& impl) : StateInputConstraint(ConstraintOrder::Linear), impl_(*impl.clone()) {}
  FrictionConeIneqConstraint(const std::unique_ptr<FrictionConeConstraint>& impl) : StateInputConstraint(ConstraintOrder::Linear), 
                                                                                    impl_(impl->clone()) {}

  ~FrictionConeIneqConstraint() override = default;
  FrictionConeIneqConstraint* clone() const override { return new FrictionConeIneqConstraint(impl_); }

  bool isActive(scalar_t time) const override;
  size_t getNumConstraints(scalar_t time) const override { return 1; };
  vector_t getValue(scalar_t time, const vector_t& state, const vector_t& input, const PreComputation& preComp) const override;
  VectorFunctionLinearApproximation getLinearApproximation(scalar_t time, const vector_t& state, const vector_t& input,
                                                           const PreComputation& preComp) const override;

  /** Sets the estimated terrain normal expressed in the world frame. */
  void setSurfaceNormalInWorld(const vector3_t& surfaceNormalInWorld);

 private:
  std::unique_ptr<FrictionConeConstraint> impl_;
};

}  // namespace legged_robot
}  // namespace ocs2
