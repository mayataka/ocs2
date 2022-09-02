#include "ocs2_sto/constraint/StoConstraintCollection.h"

namespace ocs2 {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
StoConstraintCollection::StoConstraintCollection(const StoConstraintCollection& other) : Collection<StoConstraint>(other) {}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
StoConstraintCollection* StoConstraintCollection::clone() const {
  return new StoConstraintCollection(*this);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
size_t StoConstraintCollection::getNumConstraints(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, 
                                                  const ModeSchedule& stoModeSchedule) const {
  size_t numConstraints = 0;

  // accumulate number of constraints for each constraintTerm
  for (const auto& constraintTerm : this->terms_) {
    if (constraintTerm->isActive(initTime, finalTime, referenceModeSchedule, stoModeSchedule)) {
      numConstraints += constraintTerm->getNumConstraints(initTime, finalTime, referenceModeSchedule, stoModeSchedule);
    }
  }

  return numConstraints;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t StoConstraintCollection::getValue(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, 
                                           const ModeSchedule& stoModeSchedule, const PreComputation& preComp) const {
  vector_t constraintValues;
  constraintValues.resize(getNumConstraints(initTime, finalTime, referenceModeSchedule, stoModeSchedule));

  // append vectors of constraint values from each constraintTerm
  size_t i = 0;
  for (const auto& constraintTerm : this->terms_) {
    if (constraintTerm->isActive(initTime, finalTime, referenceModeSchedule, stoModeSchedule)) {
      const auto constraintTermValues = constraintTerm->getValue(initTime, finalTime, referenceModeSchedule, stoModeSchedule, preComp);
      constraintValues.segment(i, constraintTermValues.rows()) = constraintTermValues;
      i += constraintTermValues.rows();
    }
  }

  return constraintValues;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximation StoConstraintCollection::getLinearApproximation(scalar_t initTime, scalar_t finalTime, 
                                                                                  const ModeSchedule& referenceModeSchedule,
                                                                                  const ModeSchedule& stoModeSchedule, 
                                                                                  const PreComputation& preComp) const {
  const auto numSwitchingTimes = getNumValidSwitchingTimes(initTime, finalTime, referenceModeSchedule);
  VectorFunctionLinearApproximation linearApproximation(getNumConstraints(initTime, finalTime, referenceModeSchedule, stoModeSchedule), 
                                                        numSwitchingTimes, 0);
  // append linearApproximation of each constraintTerm
  size_t i = 0;
  for (const auto& constraintTerm : this->terms_) {
    if (constraintTerm->isActive(initTime, finalTime, referenceModeSchedule, stoModeSchedule)) {
      const auto constraintTermApproximation = constraintTerm->getLinearApproximation(initTime, finalTime, referenceModeSchedule, 
                                                                                      stoModeSchedule, preComp);
      const size_t nc = constraintTermApproximation.f.rows();
      linearApproximation.f.segment(i, nc) = constraintTermApproximation.f;
      linearApproximation.dfdx.middleRows(i, nc) = constraintTermApproximation.dfdx;
      i += nc;
    }
  }

  return linearApproximation;
}

}  // namespace ocs2
