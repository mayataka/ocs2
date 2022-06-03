#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/oc_problem/OptimalControlProblem.h>

namespace ocs2 {
namespace ipm {

struct Metrics {
  using value_t = std::pair<vector_t, scalar_t>;

  // Cost
  scalar_t cost;

  // (Discretized) state equation  
  vector_t stateEquation;

  // Inequality constraints
  scalar_t slackBarrier;
  vector_t stateIneqConstraint;
  vector_t stateInputIneqConstraint;

  // Equality constraints
  vector_t stateEqConstraint;
  vector_t stateInputEqConstraint;

  // Lagrangians
  //  std::vector<value_t> stateEqLagrangian;
  //  std::vector<value_t> stateIneqLagrangian;
  //  std::vector<value_t> stateInputEqLagrangian;
  //  std::vector<value_t> stateInputIneqLagrangian;
  scalar_t stateEqLagrangian;
  scalar_t stateIneqLagrangian;
  scalar_t stateInputEqLagrangian;
  scalar_t stateInputIneqLagrangian;
};

struct MetricsCollection {
  Metrics final;
  std::vector<Metrics> preJumps;
  std::vector<Metrics> intermediates;
};

/** Exchanges the given values of Metrics */
void swap(Metrics& lhs, Metrics& rhs);

/** Clears the value of the given Metrics */
void clear(Metrics& m);

/** Exchanges the given values of MetricsCollection */
void swap(MetricsCollection& lhs, MetricsCollection& rhs);

/** Clears the value of the given MetricsCollection */
void clear(MetricsCollection& m);

}  // namespace ipm
}  // namespace ocs2
