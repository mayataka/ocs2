#include "ocs2_ipm/oc_data/Metrics.h"

namespace ocs2 {
namespace ipm {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void swap(Metrics& lhs, Metrics& rhs) {
  // Cost
  std::swap(lhs.cost, rhs.cost);

  // Inequality constraints
  lhs.stateIneqConstraint.swap(rhs.stateIneqConstraint);
  lhs.stateInputIneqConstraint.swap(rhs.stateInputIneqConstraint);

  // Equality constraints
  lhs.stateEqConstraint.swap(rhs.stateEqConstraint);
  lhs.stateInputEqConstraint.swap(rhs.stateInputEqConstraint);

  // Lagrangians
  //  lhs.stateEqLagrangian.swap(rhs.stateEqLagrangian);
  //  lhs.stateIneqLagrangian.swap(rhs.stateIneqLagrangian);
  //  lhs.stateInputEqLagrangian.swap(rhs.stateInputEqLagrangian);
  //  lhs.stateInputIneqLagrangian.swap(rhs.stateInputIneqLagrangian);

  std::swap(lhs.stateEqLagrangian, rhs.stateEqLagrangian);
  std::swap(lhs.stateIneqLagrangian, rhs.stateIneqLagrangian);
  std::swap(lhs.stateInputEqLagrangian, rhs.stateInputEqLagrangian);
  std::swap(lhs.stateInputIneqLagrangian, rhs.stateInputIneqLagrangian);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void clear(Metrics& m) {
  // Cost
  m.cost = 0.0;

  // Inequality constraints
  m.stateIneqConstraint = vector_t();
  m.stateInputIneqConstraint = vector_t();

  // Equality constraints
  m.stateEqConstraint = vector_t();
  m.stateInputEqConstraint = vector_t();

  // Lagrangians
  //  m.stateEqLagrangian.clear();
  //  m.stateIneqLagrangian.clear();
  //  m.stateInputEqLagrangian.clear();
  //  m.stateInputIneqLagrangian.clear();

  m.stateEqLagrangian = 0.0;
  m.stateIneqLagrangian = 0.0;
  m.stateInputEqLagrangian = 0.0;
  m.stateInputIneqLagrangian = 0.0;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void swap(MetricsCollection& lhs, MetricsCollection& rhs) {
  swap(lhs.final, rhs.final);
  lhs.preJumps.swap(rhs.preJumps);
  lhs.intermediates.swap(rhs.intermediates);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void clear(MetricsCollection& m) {
  clear(m.final);
  m.preJumps.clear();
  m.intermediates.clear();
}

}  // namespace ipm
}  // namespace ocs2
