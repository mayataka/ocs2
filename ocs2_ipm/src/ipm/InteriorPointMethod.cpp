#include <iostream>
#include <cassert>

#include <ocs2_ipm/ipm/InteriorPointMethod.h>

namespace ocs2 {

void initSlackDualData(const VectorFunctionLinearApproximation& constraint,
                       SlackDualData& slackDual, scalar_t barrier) {
  const auto nc = constraint.f.size();
  setBarrier(nc, barrier, slackDual);
  slackDual.slack = - constraint.f;
  for (size_t i=0; i<nc; ++i) {
    slackDual.slack[i] = std::max(slackDual.slack[i], std::sqrt(barrier));
  }
  slackDual.dual.array() = barrier / slackDual.slack.array();
}


void initInteriorPointMethodData(const VectorFunctionLinearApproximation& constraint,
                                 InteriorPointMethodData& data) {
  const auto nc = constraint.f.size();
  const auto nx = constraint.dfdx.cols();
  const auto nu = constraint.dfdu.cols();
  data.setZero(nc, nx, nu);
}


void evalPerturbedResidual(const VectorFunctionLinearApproximation& constraint,
                           const SlackDualData& slackDual, 
                           InteriorPointMethodData& data) {
  data.primalResidual = constraint.f + slackDual.slack;
  data.complementary.array() = slackDual.slack.array() * slackDual.dual.array() - slackDual.barrier;
}


void evalLagrangianDerivativesStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateConstraint, 
    const SlackDualData& slackDual, ScalarFunctionQuadraticApproximation& cost) {
  cost.dfdx.noalias() += stateConstraint.dfdx.transpose() * slackDual.dual;
}


void evalLagrangianDerivativesInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputConstraint, 
    const SlackDualData& slackDual, ScalarFunctionQuadraticApproximation& cost) {
  cost.dfdu.noalias() += inputConstraint.dfdu.transpose() * slackDual.dual;
}


void evalLagrangianDerivativesStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputConstraint, 
    const SlackDualData& slackDual, ScalarFunctionQuadraticApproximation& cost) {
  evalLagrangianDerivativesStateIneqConstraint(stateInputConstraint, slackDual, cost);
  evalLagrangianDerivativesInputIneqConstraint(stateInputConstraint, slackDual, cost);
}


void computeCondensingCoefficient(const SlackDualData& slackDual, InteriorPointMethodData& data) {
  data.cond.array() = (slackDual.dual.array()*data.primalResidual.array()-data.complementary.array()) 
                        / slackDual.slack.array();
  data.dualDivSlack.array() = slackDual.dual.array() / slackDual.slack.array();
}


void condensingStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateIneqConstraint,
    const SlackDualData& slackDual, InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  computeCondensingCoefficient(slackDual, data);
  cost.dfdx.noalias() += stateIneqConstraint.dfdx.transpose() * data.cond;
  data.linearApproximation.dfdx.noalias() = data.dualDivSlack.asDiagonal() * stateIneqConstraint.dfdx;
  cost.dfdxx.noalias() += stateIneqConstraint.dfdx.transpose() * data.linearApproximation.dfdx;
}


void condensingInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputIneqConstraint,
    const SlackDualData& slackDual, InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  computeCondensingCoefficient(slackDual, data);
  cost.dfdu.noalias() += inputIneqConstraint.dfdu.transpose() * data.cond;
  data.linearApproximation.dfdu.noalias() = data.dualDivSlack.asDiagonal() * inputIneqConstraint.dfdu;
  cost.dfduu.noalias() += inputIneqConstraint.dfdu.transpose() * data.linearApproximation.dfdu;
}


void condensingStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputIneqConstraint,
    const SlackDualData& slackDual, InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost) {
  computeCondensingCoefficient(slackDual, data);
  cost.dfdx.noalias() += stateInputIneqConstraint.dfdx.transpose() * data.cond;
  cost.dfdu.noalias() += stateInputIneqConstraint.dfdu.transpose() * data.cond;
  data.linearApproximation.dfdx.noalias() = data.dualDivSlack.asDiagonal() * stateInputIneqConstraint.dfdx;
  data.linearApproximation.dfdu.noalias() = data.dualDivSlack.asDiagonal() * stateInputIneqConstraint.dfdu;
  cost.dfdxx.noalias() += stateInputIneqConstraint.dfdx.transpose() * data.linearApproximation.dfdx;
  cost.dfdux.noalias() += stateInputIneqConstraint.dfdu.transpose() * data.linearApproximation.dfdx;
  cost.dfduu.noalias() += stateInputIneqConstraint.dfdu.transpose() * data.linearApproximation.dfdu;
}


void computeDualDirection(const SlackDualData& slackDual, InteriorPointMethodData& data) {
  data.dualDirection.array() = - (slackDual.dual.array()*data.slackDirection.array()+data.complementary.array())
                                  / slackDual.slack.array();
}


void expansionStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateIneqConstraint,
    const SlackDualData& slackDual, const vector_t& dx, InteriorPointMethodData& data) {
  data.slackDirection = - data.primalResidual;
  data.slackDirection.noalias() -= stateIneqConstraint.dfdx * dx;
  computeDualDirection(slackDual, data);
}


void expansionInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputIneqConstraint,
    const SlackDualData& slackDual, const vector_t& du, InteriorPointMethodData& data) {
  data.slackDirection = - data.primalResidual;
  data.slackDirection.noalias() -= inputIneqConstraint.dfdu * du;
  computeDualDirection(slackDual, data);
}


void expansionStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputIneqConstraint, 
    const SlackDualData& slackDual, const vector_t& dx, const vector_t& du, 
    InteriorPointMethodData& data) {
  data.slackDirection = - data.primalResidual;
  data.slackDirection.noalias() -= stateInputIneqConstraint.dfdx * dx;
  data.slackDirection.noalias() -= stateInputIneqConstraint.dfdu * du;
  computeDualDirection(slackDual, data);
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


scalar_t fractionToBoundaryPrimalStepSize(const SlackDualData& slackDual, 
                                          const InteriorPointMethodData& data,
                                          scalar_t marginRate) {
  return fractionToBoundaryStepSize(data.dim, slackDual.slack, data.slackDirection, marginRate);
}


scalar_t fractionToBoundaryDualStepSize(const SlackDualData& slackDual, 
                                        const InteriorPointMethodData& data,
                                        scalar_t marginRate) {
  return fractionToBoundaryStepSize(data.dim, slackDual.dual, data.dualDirection, marginRate);
}

}  // namespace ocs2