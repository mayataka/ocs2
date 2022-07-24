#pragma once

#include <ocs2_core/Types.h>

namespace ocs2 {

/**
 * The model data for the switching time optimization problems.
 */
struct StoModelData {
  int dim = 0;

  // Cost
  ScalarFunctionQuadraticApproximation stoCost;

  // Inequality constraints
  VectorFunctionLinearApproximation stoConstraint;
};

}  // namespace ocs2
