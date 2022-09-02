#pragma once

#include <type_traits>

#include <ocs2_core/PreComputation.h>
#include <ocs2_core/Types.h>
#include <ocs2_core/reference/ModeSchedule.h>

namespace ocs2 {

/** Sto constraint function base class */
class StoConstraint {
 public:
  StoConstraint() = default;
  virtual ~StoConstraint() = default;
  virtual StoConstraint* clone() const = 0;

  /** Check if constraint term is active */
  virtual bool isActive(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, 
                        const ModeSchedule& stoModeSchedule) const { return true; }

  /** Get the size of the constraint vector at given times */
  virtual size_t getNumConstraints(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, 
                                   const ModeSchedule& stoModeSchedule) const = 0;

  /** Get the constraint vector value */
  virtual vector_t getValue(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, 
                            const ModeSchedule& stoModeSchedule, const PreComputation& preComp) const = 0;

  /** Get the constraint linear approximation */
  virtual VectorFunctionLinearApproximation getLinearApproximation(scalar_t initTime, scalar_t finalTime, 
                                                                   const ModeSchedule& referenceModeSchedule, 
                                                                   const ModeSchedule& stoModeSchedule,
                                                                   const PreComputation& preComp) const = 0;

 protected:
  StoConstraint(const StoConstraint& rhs) = default;

 private:
};

// Template for conditional compilation using SFINAE
template <typename T>
using EnableIfStoConstraint_t = typename std::enable_if<std::is_same<T, StoConstraint>::value, bool>::type;

} // namespace ocs2