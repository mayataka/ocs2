#include <ocs2_sto/constraint/StoConstraintCollection.h>

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
size_t StoConstraintCollection::getNumConstraints(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule, 
                                                  const ModeSchedule& referenceModeSchedule) const {
  size_t numConstraints = 0;

  // accumulate number of constraints for each constraintTerm
  for (const auto& constraintTerm : this->terms_) {
    if (constraintTerm->isActive(initTime, finalTime, stoModeSchedule, referenceModeSchedule)) {
      numConstraints += constraintTerm->getNumConstraints(initTime, finalTime, stoModeSchedule, referenceModeSchedule);
    }
  }

  return numConstraints;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t StoConstraintCollection::getValue(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule, 
                                           const ModeSchedule& referenceModeSchedule, const PreComputation& preComp) const {
  vector_t constraintValues;
  constraintValues.resize(getNumConstraints(initTime, finalTime, stoModeSchedule, referenceModeSchedule));

  // append vectors of constraint values from each constraintTerm
  size_t i = 0;
  for (const auto& constraintTerm : this->terms_) {
    if (constraintTerm->isActive(initTime, finalTime, stoModeSchedule, referenceModeSchedule)) {
      const auto constraintTermValues = constraintTerm->getValue(initTime, finalTime, stoModeSchedule, referenceModeSchedule, preComp);
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
                                                                                  const ModeSchedule& stoModeSchedule, 
                                                                                  const ModeSchedule& referenceModeSchedule,
                                                                                  const PreComputation& preComp) const {
  const auto numSwitchingTimes = extractValidSwitchingTimes(initTime, finalTime, referenceModeSchedule).size();
  VectorFunctionLinearApproximation linearApproximation(getNumConstraints(initTime, finalTime, stoModeSchedule, referenceModeSchedule), 
                                                        numSwitchingTimes, 0);
  // append linearApproximation of each constraintTerm
  size_t i = 0;
  for (const auto& constraintTerm : this->terms_) {
    if (constraintTerm->isActive(initTime, finalTime, stoModeSchedule, referenceModeSchedule)) {
      const auto constraintTermApproximation = constraintTerm->getLinearApproximation(initTime, finalTime, stoModeSchedule, 
                                                                                      referenceModeSchedule, preComp);
      const size_t nc = constraintTermApproximation.f.rows();
      linearApproximation.f.segment(i, nc) = constraintTermApproximation.f;
      linearApproximation.dfdx.middleRows(i, nc) = constraintTermApproximation.dfdx;
      i += nc;
    }
  }

  return linearApproximation;
}

}  // namespace ocs2
