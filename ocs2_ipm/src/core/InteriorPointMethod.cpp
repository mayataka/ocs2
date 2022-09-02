#include "ocs2_ipm/core/InteriorPointMethod.h"

#include <iostream>
#include <cassert>

namespace ocs2 {
namespace ipm {

void initSlackDual(const VectorFunctionLinearApproximation& ineqConstraint, SlackDual& slackDual, scalar_t barrierParam) {
  assert(barrierParam > 0.0);
  const auto nc = ineqConstraint.f.size();
  setBarrier(nc, barrierParam, slackDual);
  slackDual.slack = - ineqConstraint.f;
  for (size_t i=0; i<nc; ++i) {
    slackDual.slack[i] = std::max(slackDual.slack[i], std::sqrt(barrierParam));
  }
  slackDual.dual.array() = barrierParam / slackDual.slack.array();
}


SlackDual initSlackDual(const VectorFunctionLinearApproximation& ineqConstraint, scalar_t barrierParam) {
  SlackDual slackDual;
  initSlackDual(ineqConstraint, slackDual, barrierParam);
  return slackDual;
}


void initInteriorPointMethodData(const VectorFunctionLinearApproximation& ineqConstraint, InteriorPointMethodData& ipmData) {
  const auto nc = ineqConstraint.f.size();
  const auto nx = ineqConstraint.dfdx.cols();
  const auto nu = ineqConstraint.dfdu.cols();
  ipmData.setZero(nc, nx, nu);
}


InteriorPointMethodData initInteriorPointMethodData(const VectorFunctionLinearApproximation& ineqConstraint) {
  InteriorPointMethodData ipmData;
  initInteriorPointMethodData(ineqConstraint, ipmData);
  return ipmData;
}


void evalPerturbedResidual(const VectorFunctionLinearApproximation& ineqConstraint, const SlackDual& slackDual, 
                           InteriorPointMethodData& ipmData, ScalarFunctionQuadraticApproximation& cost, scalar_t barrierParam,
                           bool stateConstraint, bool inputConstraint) {
  assert(barrierParam > 0);
  if (ipmData.dim == 0) return;

  // primal feasibility
  ipmData.primalResidual = ineqConstraint.f + slackDual.slack;
  // complementary slackness
  ipmData.complementarySlackness.array() = slackDual.slack.array() * slackDual.dual.array() - barrierParam;
  // dual feasibility
  if (stateConstraint) {
    cost.dfdx.noalias() += ineqConstraint.dfdx.transpose() * slackDual.dual;
  }
  if (inputConstraint) {
    cost.dfdu.noalias() += ineqConstraint.dfdu.transpose() * slackDual.dual;
  }
  // cost barrier
  if (slackDual.slack.size() > 0) {
    ipmData.costBarrier = - barrierParam * slackDual.slack.array().log().sum();
  }
  else {
    ipmData.costBarrier = 0.0;
  }
}


void condenseSlackDual(const VectorFunctionLinearApproximation& ineqConstraint, const SlackDual& slackDual, 
                       InteriorPointMethodData& ipmData, ScalarFunctionQuadraticApproximation& cost,
                       bool stateConstraint, bool inputConstraint) {
  if (ipmData.dim == 0) return;

  // some coefficients for condensing
  ipmData.cond.array() = (slackDual.dual.array()*ipmData.primalResidual.array()-ipmData.complementarySlackness.array()) 
                        / slackDual.slack.array();
  ipmData.dualDivSlack.array() = slackDual.dual.array() / slackDual.slack.array();
  // condensing
  if (stateConstraint) {
    cost.dfdx.noalias() += ineqConstraint.dfdx.transpose() * ipmData.cond;
    ipmData.linearApproximation.dfdx.noalias() = ipmData.dualDivSlack.asDiagonal() * ineqConstraint.dfdx;
    cost.dfdxx.noalias() += ineqConstraint.dfdx.transpose() * ipmData.linearApproximation.dfdx;
  }
  if (inputConstraint) {
    cost.dfdu.noalias() += ineqConstraint.dfdu.transpose() * ipmData.cond;
    ipmData.linearApproximation.dfdu.noalias() = ipmData.dualDivSlack.asDiagonal() * ineqConstraint.dfdu;
    cost.dfduu.noalias() += ineqConstraint.dfdu.transpose() * ipmData.linearApproximation.dfdu;
  }
  if (stateConstraint && inputConstraint) {
    cost.dfdux.noalias() += ineqConstraint.dfdu.transpose() * ipmData.linearApproximation.dfdx;
  }
}


void expandDual(const SlackDual& slackDual, const InteriorPointMethodData& ipmData, 
                SlackDualDirection& slackDualDirection) {
  if (ipmData.dim == 0) return;

  slackDualDirection.dualDirection.array() 
      = - (slackDual.dual.array()*slackDualDirection.slackDirection.array()+ipmData.complementarySlackness.array())
            / slackDual.slack.array();
}


void expandSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                     const SlackDual& slackDual, const InteriorPointMethodData& ipmData, 
                     const vector_t& dx, const vector_t& du, SlackDualDirection& slackDualDirection) {
  if (ipmData.dim == 0) return;

  slackDualDirection.slackDirection = - ipmData.primalResidual;
  slackDualDirection.slackDirection.noalias() -= ineqConstraint.dfdx * dx;
  slackDualDirection.slackDirection.noalias() -= ineqConstraint.dfdu * du;
  expandDual(slackDual, ipmData, slackDualDirection);
}


void expandSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                     const SlackDual& slackDual, const InteriorPointMethodData& ipmData, 
                     const vector_t& dx, SlackDualDirection& slackDualDirection) {
  if (ipmData.dim == 0) return;

  slackDualDirection.slackDirection = - ipmData.primalResidual;
  slackDualDirection.slackDirection.noalias() -= ineqConstraint.dfdx * dx;
  expandDual(slackDual, ipmData, slackDualDirection);
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


scalar_t fractionToBoundaryPrimalStepSize(const SlackDual& slackDual, const SlackDualDirection& slackDualDirection,
                                          const InteriorPointMethodData& ipmData, scalar_t marginRate) {
  return fractionToBoundaryStepSize(ipmData.dim, slackDual.slack, slackDualDirection.slackDirection, marginRate);
}


scalar_t fractionToBoundaryDualStepSize(const SlackDual& slackDual, const SlackDualDirection& slackDualDirection,
                                        const InteriorPointMethodData& ipmData, scalar_t marginRate) {
  return fractionToBoundaryStepSize(ipmData.dim, slackDual.dual, slackDualDirection.dualDirection, marginRate);
}

}  // namespace ipm
}  // namespace ocs2