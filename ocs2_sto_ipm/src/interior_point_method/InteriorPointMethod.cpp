#include <iostream>

#include <ocs2_sto_ipm/interior_point_method/InteriorPointMethod.h>

namespace ocs2 {
namespace sto_ipm {

void initInteriorPointMethodData(const VectorFunctionLinearApproximation& constraint,
                                 InteriorPointMethodData& data, scalar_t barrier) {
  const auto nc = constraint.f.size();
  const auto nx = constraint.dfdx.cols();
  const auto nu = constraint.dfdu.cols();
  data.setBarrier(nc, nx, nu, barrier);
  data.slack = - constraint.f;
  for (size_t i=0; i<nc; ++i) {
    data.slack[i] = std::max(data.slack[i], std::sqrt(barrier));
  }
  data.dual.array() = barrier / data.slack.array();
}


void evalPerturbedResidual(const VectorFunctionLinearApproximation& constraint,
                           InteriorPointMethodData& data) {
  data.primalResidual = constraint.f + data.slack;
  data.complementary.array() = data.slack.array() * data.dual.array() - data.barrier;
}


void evalLagrangianDerivativesStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateConstraint, 
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  cost.dfdx.noalias() += stateConstraint.dfdx.transpose() * data.dual;
}


void evalLagrangianDerivativesInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputConstraint, 
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  cost.dfdu.noalias() += inputConstraint.dfdu.transpose() * data.dual;
}


void evalLagrangianDerivativesStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputConstraint, 
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  evalLagrangianDerivativesStateIneqConstraint(stateInputConstraint, data, cost);
  evalLagrangianDerivativesInputIneqConstraint(stateInputConstraint, data, cost);
}


void computeCondensingCoefficient(InteriorPointMethodData& data) {
  data.cond.array() = (data.dual.array()*data.primalResidual.array()-data.complementary.array()) 
                        / data.slack.array();
}


void condensingStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateIneqConstraint,
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  computeCondensingCoefficient(data);
  cost.dfdx.noalias() += stateIneqConstraint.dfdx.transpose() * data.cond;
  data.cond.array() = data.dual.array() / data.slack.array();
  data.linearApproximation.dfdx.noalias() = data.cond.asDiagonal() * stateIneqConstraint.dfdx;
  cost.dfdxx.noalias() += stateIneqConstraint.dfdx.transpose() * data.linearApproximation.dfdx;
}


void condensingInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputIneqConstraint,
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  computeCondensingCoefficient(data);
  cost.dfdu.noalias() += inputIneqConstraint.dfdu.transpose() * data.cond;
  data.cond.array() = data.dual.array() / data.slack.array();
  data.linearApproximation.dfdu.noalias() = data.cond.asDiagonal() * inputIneqConstraint.dfdu;
  cost.dfduu.noalias() += inputIneqConstraint.dfdu.transpose() * data.linearApproximation.dfdu;
}


void condensingStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputIneqConstraint,
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  computeCondensingCoefficient(data);
  cost.dfdx.noalias() += stateInputIneqConstraint.dfdx.transpose() * data.cond;
  cost.dfdu.noalias() += stateInputIneqConstraint.dfdu.transpose() * data.cond;
  data.cond.array() = data.dual.array() / data.slack.array();
  data.linearApproximation.dfdx.noalias() = data.cond.asDiagonal() * stateInputIneqConstraint.dfdx;
  data.linearApproximation.dfdu.noalias() = data.cond.asDiagonal() * stateInputIneqConstraint.dfdu;
  cost.dfdxx.noalias() += stateInputIneqConstraint.dfdx.transpose() * data.linearApproximation.dfdx;
  cost.dfdux.noalias() += stateInputIneqConstraint.dfdu.transpose() * data.linearApproximation.dfdx;
  cost.dfduu.noalias() += stateInputIneqConstraint.dfdu.transpose() * data.linearApproximation.dfdu;
}


}  // namespace sto_ipm
}  // namespace ocs2