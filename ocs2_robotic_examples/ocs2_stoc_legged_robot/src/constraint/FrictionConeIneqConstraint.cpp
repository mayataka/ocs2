#include "ocs2_stoc_legged_robot/constraint/FrictionConeIneqConstraint.h"

#include <ocs2_centroidal_model/AccessHelperFunctions.h>

namespace ocs2 {
namespace legged_robot {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
FrictionConeIneqConstraint::FrictionConeIneqConstraint(const SwitchedModelReferenceManager& referenceManager, Config config,
                                                       size_t contactPointIndex, CentroidalModelInfo info)
    : StateInputConstraint(ConstraintOrder::Linear),
      impl_(std::unique_ptr<FrictionConeConstraint>(new FrictionConeConstraint(referenceManager, config, contactPointIndex, info))) {}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void FrictionConeIneqConstraint::setSurfaceNormalInWorld(const vector3_t& surfaceNormalInWorld) {
  impl_->setSurfaceNormalInWorld(surfaceNormalInWorld);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
bool FrictionConeIneqConstraint::isActive(scalar_t time) const {
  return impl_->isActive(time);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t FrictionConeIneqConstraint::getValue(scalar_t time, const vector_t& state, const vector_t& input,
                                              const PreComputation& preComp) const {
  const vector_t value = - impl_->getValue(time, state, input, preComp);
  return value;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximation FrictionConeIneqConstraint::getLinearApproximation(scalar_t time, const vector_t& state,
                                                                                     const vector_t& input,
                                                                                     const PreComputation& preComp) const {
  VectorFunctionLinearApproximation linearApproximation = impl_->getLinearApproximation(time, state, input, preComp);
  linearApproximation.f.array() *= -1.0;
  linearApproximation.dfdx.array() *= -1.0;
  linearApproximation.dfdu.array() *= -1.0;
  return linearApproximation;
}

}  // namespace legged_robot
}  // namespace ocs2
