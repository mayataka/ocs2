#pragma once

#include <unordered_map>
#include <algorithm>

#include <ocs2_core/PreComputation.h>
#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>
#include <ocs2_core/NumericTraits.h>

#include "ocs2_sto/constraint/StoConstraint.h"
#include "ocs2_sto/ValidModeSchedule.h"


namespace ocs2 {

class MinimumDwellTimeConstraint : public StoConstraint {
 public:
  MinimumDwellTimeConstraint(const std::unordered_map<size_t, scalar_t>& minimumDwellTimesMap,
                             const scalar_t minimumDwellTime=numeric_traits::limitEpsilon<scalar_t>());
  MinimumDwellTimeConstraint() = default;
  virtual ~MinimumDwellTimeConstraint() = default;
  virtual MinimumDwellTimeConstraint* clone() const;

  /** Get the size of the constraint vector at given times */
  virtual size_t getNumConstraints(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule,  
                                   const ModeSchedule& referenceModeSchedule) const override;

  /** Get the constraint vector value */
  virtual vector_t getValue(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule, 
                            const ModeSchedule& referenceModeSchedule, const PreComputation& preComp) const override;

  /** Get the constraint linear approximation */
  virtual VectorFunctionLinearApproximation getLinearApproximation(scalar_t initTime, scalar_t finalTime, 
                                                                   const ModeSchedule& stoModeSchedule,
                                                                   const ModeSchedule& referenceModeSchedule,
                                                                   const PreComputation& preComp) const override;

 protected:
  MinimumDwellTimeConstraint(const MinimumDwellTimeConstraint& rhs) = default;

 private:
  std::unordered_map<size_t, scalar_t> minimumDwellTimesMap_;
  scalar_t minimumDwellTime_;

};

} // namespace ocs2