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
size_t StoConstraintCollection::getNumConstraints(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime, 
                                                  const ModeSchedule& modeSchedule) const {
  size_t numConstraints = 0;

  // accumulate number of constraints for each constraintTerm
  for (const auto& constraintTerm : this->terms_) {
    if (constraintTerm->isActive(initTime, switchingTimes, finalTime, modeSchedule)) {
      numConstraints += constraintTerm->getNumConstraints(initTime, switchingTimes, finalTime, modeSchedule);
    }
  }

  return numConstraints;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
vector_t StoConstraintCollection::getValue(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime, 
                                           const ModeSchedule& modeSchedule, const PreComputation& preComp) const {
  vector_t constraintValues;
  constraintValues.resize(getNumConstraints(initTime, switchingTimes, finalTime, modeSchedule));

  // append vectors of constraint values from each constraintTerm
  size_t i = 0;
  for (const auto& constraintTerm : this->terms_) {
    if (constraintTerm->isActive(initTime, switchingTimes, finalTime, modeSchedule)) {
      const auto constraintTermValues = constraintTerm->getValue(initTime, switchingTimes, finalTime, modeSchedule, preComp);
      constraintValues.segment(i, constraintTermValues.rows()) = constraintTermValues;
      i += constraintTermValues.rows();
    }
  }

  return constraintValues;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximation StoConstraintCollection::getLinearApproximation(scalar_t initTime, const vector_t& switchingTimes, 
                                                                                  scalar_t finalTime, const ModeSchedule& modeSchedule,
                                                                                  const PreComputation& preComp) const {
  VectorFunctionLinearApproximation linearApproximation(getNumConstraints(initTime, switchingTimes, finalTime, modeSchedule), 
                                                        switchingTimes.size(), 0);

  // append linearApproximation of each constraintTerm
  size_t i = 0;
  for (const auto& constraintTerm : this->terms_) {
    if (constraintTerm->isActive(initTime, switchingTimes, finalTime, modeSchedule)) {
      const auto constraintTermApproximation = constraintTerm->getLinearApproximation(initTime, switchingTimes, finalTime, modeSchedule, preComp);
      const size_t nc = constraintTermApproximation.f.rows();
      linearApproximation.f.segment(i, nc) = constraintTermApproximation.f;
      linearApproximation.dfdx.middleRows(i, nc) = constraintTermApproximation.dfdx;
      i += nc;
    }
  }

  return linearApproximation;
}

}  // namespace ocs2
