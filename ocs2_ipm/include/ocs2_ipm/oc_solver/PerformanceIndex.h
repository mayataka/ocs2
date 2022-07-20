#pragma once

#include <iomanip>
#include <ostream>

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_ipm/oc_data/IpmData.h>
#include <ocs2_oc/oc_solver/PerformanceIndex.h>

namespace ocs2 {
namespace ipm {

/**
 * Defines the performance indices for a rollout
 */
struct PerformanceIndex {
  /** The total cost. */
  scalar_t cost = 0.0;

  /** The total cost of the barrier function on the slack variables. */
  scalar_t costBarrier = 0.0;

  /** Sum of Squared Error (SSE) of system dynamics violation */
  scalar_t dynamicsViolationSSE = 0.0;

  /** Sum of Squared Error (SSE) of equality constraints:
   * - Final: squared norm of violation in state equality constraints
   * - PreJumps: sum of squared norm of violation in state equality constraints
   * - Intermediates: Integral of squared norm violation in state/state-input equality constraints
   */
  scalar_t equalityConstraintsSSE = 0.0;

  /** Sum of Squared Error (SSE) of inequality constraints:
   * - Final: squared norm of violation in state inequality constraints
   * - PreJumps: sum of squared norm of violation in state inequality constraints
   * - Intermediates: Integral of squared norm violation in state/state-input inequality constraints
   */
  scalar_t inequalityConstraintsSSE = 0.0;

  /** Sum of equality Lagrangians:
   * 
   * - Final: penalty for violation in state equality constraints
   * - PreJumps: penalty for violation in state equality constraints
   * - Intermediates: penalty for violation in state/state-input equality constraints
   */
  scalar_t equalityLagrangian = 0.0;

  /** Sum of inequality Lagrangians:
   * - Final: penalty for violation in state inequality constraints
   * - PreJumps: penalty for violation in state inequality constraints
   * - Intermediates: penalty for violation in state/state-input inequality constraints
   */
  scalar_t inequalityLagrangian = 0.0;

  /** Sum of Squared Error (SSE) of system dual dynamics violation */
  scalar_t dualDynamicsViolationSSE = 0.0;

  /** Sum of Squared Error (SSE) of system dual control violation */
  scalar_t dualControlViolationSSE = 0.0;

  /** Sum of Squared Error (SSE) of complementary slackness */
  scalar_t complementarySlacknessSSE = 0.0;

  /** Add performance indices */
  PerformanceIndex& operator+=(const PerformanceIndex& rhs) {
    this->cost += rhs.cost;
    this->costBarrier += rhs.costBarrier;
    this->dynamicsViolationSSE += rhs.dynamicsViolationSSE;
    this->equalityConstraintsSSE += rhs.equalityConstraintsSSE;
    this->inequalityConstraintsSSE += rhs.inequalityConstraintsSSE;
    this->equalityLagrangian += rhs.equalityLagrangian;
    this->inequalityLagrangian += rhs.inequalityLagrangian;
    this->dualDynamicsViolationSSE += rhs.dualDynamicsViolationSSE;
    this->dualControlViolationSSE += rhs.dualControlViolationSSE;
    this->complementarySlacknessSSE += rhs.complementarySlacknessSSE;
    return *this;
  }
};

PerformanceIndex fromModelData(const ModelData& modelData) {
  PerformanceIndex performanceIndex;
  performanceIndex.cost = modelData.cost.f;
  performanceIndex.dynamicsViolationSSE = modelData.dynamics.f.squaredNorm();
  if (modelData.stateEqConstraint.f.size() > 0)
    performanceIndex.equalityConstraintsSSE += modelData.stateEqConstraint.f.squaredNorm();
  if (modelData.stateInputEqConstraint.f.size() > 0)
    performanceIndex.equalityConstraintsSSE += modelData.stateInputEqConstraint.f.squaredNorm();
  if (modelData.cost.dfdx.size() > 0)
    performanceIndex.dualDynamicsViolationSSE = modelData.cost.dfdx.squaredNorm();
  if (modelData.cost.dfdu.size() > 0)
    performanceIndex.dualControlViolationSSE = modelData.cost.dfdu.squaredNorm();
  return performanceIndex;
}

PerformanceIndex fromIpmData(const IpmData& ipmData) {
  PerformanceIndex performanceIndex;
  if (ipmData.dataStateIneqConstraint.primalResidual.size() > 0)
    performanceIndex.inequalityConstraintsSSE += ipmData.dataStateIneqConstraint.primalResidual.squaredNorm();
  if (ipmData.dataStateInputIneqConstraint.primalResidual.size() > 0)
    performanceIndex.inequalityConstraintsSSE += ipmData.dataStateInputIneqConstraint.primalResidual.squaredNorm();
  if (ipmData.dataStateIneqConstraint.complementarySlackness.size() > 0)
    performanceIndex.complementarySlacknessSSE += ipmData.dataStateIneqConstraint.complementarySlackness.squaredNorm();
  if (ipmData.dataStateInputIneqConstraint.complementarySlackness.size() > 0)
    performanceIndex.complementarySlacknessSSE += ipmData.dataStateInputIneqConstraint.complementarySlackness.squaredNorm();
  performanceIndex.costBarrier = ipmData.dataStateIneqConstraint.costBarrier 
                                  + ipmData.dataStateInputIneqConstraint.costBarrier;
  return performanceIndex;
}

::ocs2::PerformanceIndex convert(const ::ocs2::ipm::PerformanceIndex& ipmPerformanceIndex) {
  ::ocs2::PerformanceIndex performanceIndex;
  performanceIndex.cost = ipmPerformanceIndex.cost + ipmPerformanceIndex.costBarrier;
  performanceIndex.dynamicsViolationSSE = ipmPerformanceIndex.dynamicsViolationSSE;
  performanceIndex.equalityConstraintsSSE = ipmPerformanceIndex.equalityConstraintsSSE;
                                              + ipmPerformanceIndex.inequalityConstraintsSSE;
  performanceIndex.equalityLagrangian = ipmPerformanceIndex.equalityLagrangian;
  performanceIndex.inequalityLagrangian = ipmPerformanceIndex.inequalityLagrangian;
  return performanceIndex;
}

inline PerformanceIndex operator+(PerformanceIndex lhs, const PerformanceIndex& rhs) {
  lhs += rhs;  // Copied lhs, add rhs to it.
  return lhs;
}

/** Swaps performance indices */
inline void swap(PerformanceIndex& lhs, PerformanceIndex& rhs) {
  std::swap(lhs.cost, rhs.cost);
  std::swap(lhs.costBarrier, rhs.costBarrier);
  std::swap(lhs.dynamicsViolationSSE, rhs.dynamicsViolationSSE);
  std::swap(lhs.equalityConstraintsSSE, rhs.equalityConstraintsSSE);
  std::swap(lhs.inequalityConstraintsSSE, rhs.inequalityConstraintsSSE);
  std::swap(lhs.equalityLagrangian, rhs.equalityLagrangian);
  std::swap(lhs.inequalityLagrangian, rhs.inequalityLagrangian);
  std::swap(lhs.dualDynamicsViolationSSE, rhs.dualDynamicsViolationSSE);
  std::swap(lhs.dualControlViolationSSE, rhs.dualControlViolationSSE);
  std::swap(lhs.complementarySlacknessSSE, rhs.complementarySlacknessSSE);
}

inline std::ostream& operator<<(std::ostream& stream, const PerformanceIndex& performanceIndex) {
  const size_t tabSpace = 12;
  const auto indentation = stream.width();
  stream << std::left;  // fill from left

  stream << std::setw(indentation) << "";
  stream << "Cost:                        " << std::setw(tabSpace) << performanceIndex.cost;
  stream << "Cost Barrier:                " << std::setw(tabSpace) << performanceIndex.costBarrier << '\n';

  stream << std::setw(indentation) << "";
  stream << "Dynamics violation SSE:      " << std::setw(tabSpace) << performanceIndex.dynamicsViolationSSE;
  stream << "Equality constraints SSE:    " << std::setw(tabSpace) << performanceIndex.equalityConstraintsSSE;
  stream << "Inequality constraints SSE:  " << std::setw(tabSpace) << performanceIndex.inequalityConstraintsSSE << '\n';

  stream << std::setw(indentation) << "";
  stream << "Equality Lagrangian:         " << std::setw(tabSpace) << performanceIndex.equalityLagrangian;
  stream << "Inequality Lagrangian:       " << std::setw(tabSpace) << performanceIndex.inequalityLagrangian << '\n';

  stream << std::setw(indentation) << "";
  stream << "Dual dynamics violation SSE: " << std::setw(tabSpace) << performanceIndex.dualDynamicsViolationSSE;
  stream << "Dual control violation SSE:  " << std::setw(tabSpace) << performanceIndex.dualControlViolationSSE;
  stream << "Complementary Slackness SSE: " << std::setw(tabSpace) << performanceIndex.complementarySlacknessSSE;

  return stream;
}

}  // namespace ipm
}  // namespace ocs2
