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
  data.dualDivSlack.array() = data.dual.array() / data.slack.array();
}


void condensingStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateIneqConstraint,
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  computeCondensingCoefficient(data);
  cost.dfdx.noalias() += stateIneqConstraint.dfdx.transpose() * data.cond;
  data.linearApproximation.dfdx.noalias() = data.dualDivSlack.asDiagonal() * stateIneqConstraint.dfdx;
  cost.dfdxx.noalias() += stateIneqConstraint.dfdx.transpose() * data.linearApproximation.dfdx;
}


void condensingInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputIneqConstraint,
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  computeCondensingCoefficient(data);
  cost.dfdu.noalias() += inputIneqConstraint.dfdu.transpose() * data.cond;
  data.linearApproximation.dfdu.noalias() = data.dualDivSlack.asDiagonal() * inputIneqConstraint.dfdu;
  cost.dfduu.noalias() += inputIneqConstraint.dfdu.transpose() * data.linearApproximation.dfdu;
}


void condensingStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputIneqConstraint,
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  computeCondensingCoefficient(data);
  cost.dfdx.noalias() += stateInputIneqConstraint.dfdx.transpose() * data.cond;
  cost.dfdu.noalias() += stateInputIneqConstraint.dfdu.transpose() * data.cond;
  data.linearApproximation.dfdx.noalias() = data.dualDivSlack.asDiagonal() * stateInputIneqConstraint.dfdx;
  data.linearApproximation.dfdu.noalias() = data.dualDivSlack.asDiagonal() * stateInputIneqConstraint.dfdu;
  cost.dfdxx.noalias() += stateInputIneqConstraint.dfdx.transpose() * data.linearApproximation.dfdx;
  cost.dfdux.noalias() += stateInputIneqConstraint.dfdu.transpose() * data.linearApproximation.dfdx;
  cost.dfduu.noalias() += stateInputIneqConstraint.dfdu.transpose() * data.linearApproximation.dfdu;
}


void computeDualDirection(InteriorPointMethodData& data) {
  data.dualDirection.array() = - (data.dual.array()*data.slackDirection.array()+data.complementary.array())
                                  / data.slack.array();
}


void expansionStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateIneqConstraint,
    const vector_t& dx, InteriorPointMethodData& data) {
  data.slackDirection = - data.primalResidual;
  data.slackDirection.noalias() -= stateIneqConstraint.dfdx * dx;
  computeDualDirection(data);
}


void expansionInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputIneqConstraint,
    const vector_t& du, InteriorPointMethodData& data) {
  data.slackDirection = - data.primalResidual;
  data.slackDirection.noalias() -= inputIneqConstraint.dfdu * du;
  computeDualDirection(data);
}


void expansionStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputIneqConstraint, 
    const vector_t& dx, const vector_t& du, InteriorPointMethodData& data) {
  data.slackDirection = - data.primalResidual;
  data.slackDirection.noalias() -= stateInputIneqConstraint.dfdx * dx;
  data.slackDirection.noalias() -= stateInputIneqConstraint.dfdu * du;
  computeDualDirection(data);
}


scalar_t fractionToBoundaryStepSize(size_t dim, const vector_t& v, const vector_t& dv,
                                    scalar_t marginRate) {
  assert(marginRate > 0);
  assert(marginRate <= 1);
  scalar_t minFractionToBoundary = 1.0;
  for (size_t i=0; i<dim; ++i) {
    const scalar_t fractionToBoundary
        = - marginRate * (v.coeff(i)/dv.coeff(i));
    if (fractionToBoundary > 0.0 && fractionToBoundary < 1.0) {
      if (fractionToBoundary < minFractionToBoundary) {
        minFractionToBoundary = fractionToBoundary;
      }
    }
  }
  assert(minFractionToBoundary > 0);
  assert(minFractionToBoundary <= 1);
  return minFractionToBoundary;
}


scalar_t fractionToBoundaryPrimalStepSize(const InteriorPointMethodData& data,
                                          scalar_t marginRate) {
  return fractionToBoundaryStepSize(data.dim, data.slack, data.slackDirection, marginRate);
}


scalar_t fractionToBoundaryDualStepSize(const InteriorPointMethodData& data,
                                        scalar_t marginRate) {
  return fractionToBoundaryStepSize(data.dim, data.dual, data.dualDirection, marginRate);
}

}  // namespace sto_ipm
}  // namespace ocs2