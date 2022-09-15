#pragma once

#include <unordered_map>
#include <algorithm>

#include <ocs2_core/PreComputation.h>
#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>
#include <ocs2_core/NumericTraits.h>

#include "ocs2_sto/constraint/MinimumDwellTimeConstraint.h"


namespace ocs2 {

/** Minimum dwell-time constraint term for MPC */
class MPC_MinimumDwellTimeConstraint : public MinimumDwellTimeConstraint {
 public:
  /**
   * Constructor for the minimum dwell-time constraint for MPC.
   * @param [in] minimumDwellTimeAtInitialPhase: Minimum dwell time at the initial phase.
   * @param [in] minimumDwellTimesMap: Minimum dwell times for modes.
   * @param [in] minimumDwellTime: Minimum dwell time for the modes that are not specified by minimumDwellTimesMap. 
   */
  MPC_MinimumDwellTimeConstraint(const scalar_t minimumDwellTimeAtInitialPhase,
                                 const std::unordered_map<size_t, scalar_t>& minimumDwellTimesMap,
                                 const scalar_t minimumDwellTime=numeric_traits::limitEpsilon<scalar_t>());
  MPC_MinimumDwellTimeConstraint() = default;
  virtual ~MPC_MinimumDwellTimeConstraint() = default;
  virtual MPC_MinimumDwellTimeConstraint* clone() const;

 protected:
  MPC_MinimumDwellTimeConstraint(const MPC_MinimumDwellTimeConstraint& rhs) = default;

  virtual vector_t getMinimumDwellTimes(const ModeSchedule& validModeSchedule) const override;

 private:
  std::unordered_map<size_t, scalar_t> minimumDwellTimesMap_;
  scalar_t minimumDwellTime_, minimumDwellTimeAtInitialPhase_;
};

} // namespace ocs2