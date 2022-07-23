#pragma once

#include <type_traits>

#include <ocs2_core/PreComputation.h>
#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>

namespace ocs2 {

class StoCost {
 public:
  StoCost() = default;
  virtual ~StoCost() = default;
  virtual StoCost* clone() const = 0;

  /** Check if cost term is active */
  virtual bool isActive(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule, 
                        const ModeSchedule& referenceModeSchedule) const { return true; }

  /** Get cost term value */
  virtual scalar_t getValue(scalar_t initTime, scalar_t finalTime, const ModeSchedule& stoModeSchedule, 
                            const ModeSchedule& referenceModeSchedule, const PreComputation& preComp) const = 0;

  /** Get cost term quadratic approximation */
  virtual ScalarFunctionQuadraticApproximation getQuadraticApproximation(scalar_t initTime, scalar_t finalTime,  
                                                                         const ModeSchedule& stoModeSchedule, 
                                                                         const ModeSchedule& referenceModeSchedule, 
                                                                         const PreComputation& preComp) const = 0;

 protected:
  StoCost(const StoCost& rhs) = default;
};

// Template for conditional compilation using SFINAE
template <typename T>
using EnableIfStoCost_t = typename std::enable_if<std::is_same<T, StoCost>::value, bool>::type;

}  // namespace ocs2
