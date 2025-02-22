#pragma once

#include "ocs2_sto/cost/StoCost.h"
#include "ocs2_sto/ValidModeSchedule.h"

namespace ocs2 {

/** Quadratic Sto cost term */
class QuadraticStoCost : public StoCost {
 public:
  /**
   * Constructor for the quadratic cost function defined as the following:
   * \f$ \l = 0.5(x-x_{n})' (Q*I_{n}) (x-x_{n}) \f$. (x: a swithcing times vector)
   * @param [in] Q: \f$ Q \f$
   */
  explicit QuadraticStoCost(scalar_t Q);
  ~QuadraticStoCost() override = default;
  QuadraticStoCost* clone() const override;

  /** Get cost term value */
  scalar_t getValue(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, 
                    const ModeSchedule& stoModeSchedule, const PreComputation& preComp) const override;

  /** Get cost term quadratic approximation */
  ScalarFunctionQuadraticApproximation getQuadraticApproximation(scalar_t initTime, scalar_t finalTime, 
                                                                 const ModeSchedule& referenceModeSchedule,
                                                                 const ModeSchedule& stoModeSchedule, 
                                                                 const PreComputation& preComp) const override;

 protected:
  QuadraticStoCost(const QuadraticStoCost& rhs) = default;

  /** Computes the state deviation for the nominal switching times.
   * This method can be overwritten if desired mode schedule has a different dimensions. */
  virtual vector_t getSwitchingTimeDeviation(scalar_t initTime, scalar_t finalTime, const ModeSchedule& referenceModeSchedule, 
                                             const ModeSchedule& stoModeSchedule) const;

 private:
  scalar_t Q_;
};

}  // namespace ocs2
