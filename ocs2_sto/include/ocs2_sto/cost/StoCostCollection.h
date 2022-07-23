#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_core/misc/Collection.h>
#include <ocs2_core/reference/ModeSchedule.h>

#include "ocs2_sto/cost/StoCost.h"

namespace ocs2 {

/**
 * Cost function combining a collection of cost terms.
 *
 * This class collects a variable number of cost terms and provides methods to get the
 * summed cost values and quadratic approximations. Each cost term can be accessed through its
 * string name and can be activated or deactivated.
 */
class StoCostCollection : public Collection<StoCost> {
 public:
  StoCostCollection() = default;
  virtual ~StoCostCollection() = default;
  virtual StoCostCollection* clone() const;

  /** Get state-only cost value */
  virtual scalar_t getValue(scalar_t initTime, const vector_t& switchingTimes, scalar_t finalTime, const ModeSchedule& modeSchedule, 
                            const PreComputation& preComp) const;

  /** Get state-only cost quadratic approximation */
  virtual ScalarFunctionQuadraticApproximation getQuadraticApproximation(scalar_t initTime, const vector_t& switchingTimes, 
                                                                         scalar_t finalTime, const ModeSchedule& modeSchedule, 
                                                                         const PreComputation& preComp) const;

 protected:
  /** Copy constructor */
  StoCostCollection(const StoCostCollection& other);
};

}  // namespace ocs2
