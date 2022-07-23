#pragma once

#include <ocs2_core/PreComputation.h>
#include <ocs2_core/Types.h>
#include <ocs2_core/misc/Collection.h>

#include "ocs2_sto/constraint/StoConstraint.h"

namespace ocs2 {

/**
 * Constraint collection class
 *
 * This class collects a variable number of constraint functions and provides methods to get the
 * concatenated constraint vectors and approximations. Each constraint can be accessed through its
 * string name and can be activated or deactivated.
 */
class StoConstraintCollection : public Collection<StoConstraint> {
 public:
  StoConstraintCollection() = default;
  ~StoConstraintCollection() override = default;
  StoConstraintCollection* clone() const override;

  /** Get the size of the constraint vector at given times */
  virtual size_t getNumConstraints(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime, 
                                   const ModeSchedule& modeSchedule) const;

  /** Get the constraint vector value */
  virtual vector_t getValue(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime, const ModeSchedule& modeSchedule, 
                            const PreComputation& preComp) const;

  /** Get the constraint linear approximation */
  virtual VectorFunctionLinearApproximation getLinearApproximation(scalar_t initTime, const vector_t& switchingTimes, 
                                                                   scalar_t finalTime, const ModeSchedule& modeSchedule,
                                                                   const PreComputation& preComp) const;

 protected:
  /** Copy constructor */
  StoConstraintCollection(const StoConstraintCollection& other);
};

}  // namespace ocs2