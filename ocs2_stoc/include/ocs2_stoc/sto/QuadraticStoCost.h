#pragma once

#include "ocs2_stoc/sto/StoCost.h"

namespace ocs2 {

/** Quadratic state-only cost term */
class QuadraticStoCost : public StoCost {
 public:
  /**
   * Constructor for the quadratic cost function defined as the following:
   * \f$ \l = 0.5(x-x_{n})' Q (x-x_{n}) \f$. (x: a swithcing times vector)
   * @param [in] Q: \f$ Q \f$
   */
  explicit QuadraticStoCost(matrix_t Q);
  ~QuadraticStoCost() override = default;
  QuadraticStoCost* clone() const override;

  /** Get cost term value */
  scalar_t getValue(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime, const ModeSchedule& modeSchedule, 
                    const PreComputation&) const override;

  /** Get cost term quadratic approximation */
  ScalarFunctionQuadraticApproximation getQuadraticApproximation(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime,  
                                                                 const ModeSchedule& modeSchedule, const PreComputation&) const override;

 protected:
  QuadraticStoCost(const QuadraticStoCost& rhs) = default;

  /** Computes the state deviation for the nominal switching times.
   * This method can be overwritten if desired mode schedule has a different dimensions. */
  virtual vector_t getSwitchingTimeDeviation(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime,  
                                             const ModeSchedule& modeSchedule) const;

 private:
  matrix_t Q_;
};

}  // namespace ocs2
