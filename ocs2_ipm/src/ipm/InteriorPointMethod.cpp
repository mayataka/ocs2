#include <iostream>
#include <cassert>

#include <ocs2_ipm/ipm/InteriorPointMethod.h>

namespace ocs2 {

void initSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                   SlackDual& slackDual, scalar_t barrier) {
  const auto nc = ineqConstraint.f.size();
  setBarrier(nc, barrier, slackDual);
  slackDual.slack = - ineqConstraint.f;
  for (size_t i=0; i<nc; ++i) {
    slackDual.slack[i] = std::max(slackDual.slack[i], std::sqrt(barrier));
  }
  slackDual.dual.array() = barrier / slackDual.slack.array();
}


void initInteriorPointMethodData(const VectorFunctionLinearApproximation& ineqConstraint,
                                 InteriorPointMethodData& data) {
  const auto nc = ineqConstraint.f.size();
  const auto nx = ineqConstraint.dfdx.cols();
  const auto nu = ineqConstraint.dfdu.cols();
  data.setZero(nc, nx, nu);
}


void evalPerturbedResidual(const VectorFunctionLinearApproximation& ineqConstraint,
                           const SlackDual& slackDual, 
                           InteriorPointMethodData& data,
                           ScalarFunctionQuadraticApproximation& cost,
                           bool stateConstraint, bool inputConstraint) {
  // primal feasibility
  data.primalResidual = ineqConstraint.f + slackDual.slack;
  // complementary
  data.complementary.array() = slackDual.slack.array() * slackDual.dual.array() - slackDual.barrier;
  // dual feasibility
  if (stateConstraint) {
    cost.dfdx.noalias() += ineqConstraint.dfdx.transpose() * slackDual.dual;
  }
  if (inputConstraint) {
    cost.dfdu.noalias() += ineqConstraint.dfdu.transpose() * slackDual.dual;
  }
}


void condenseSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                       const SlackDual& slackDual, InteriorPointMethodData& data, 
                       ScalarFunctionQuadraticApproximation& cost,
                       bool stateConstraint, bool inputConstraint) {
  // some coefficients for condensing
  data.cond.array() = (slackDual.dual.array()*data.primalResidual.array()-data.complementary.array()) 
                        / slackDual.slack.array();
  data.dualDivSlack.array() = slackDual.dual.array() / slackDual.slack.array();
  // condensing
  if (stateConstraint) {
    cost.dfdx.noalias() += ineqConstraint.dfdx.transpose() * data.cond;
    data.linearApproximation.dfdx.noalias() = data.dualDivSlack.asDiagonal() * ineqConstraint.dfdx;
    cost.dfdxx.noalias() += ineqConstraint.dfdx.transpose() * data.linearApproximation.dfdx;
  }
  if (inputConstraint) {
    cost.dfdu.noalias() += ineqConstraint.dfdu.transpose() * data.cond;
    data.linearApproximation.dfdu.noalias() = data.dualDivSlack.asDiagonal() * ineqConstraint.dfdu;
    cost.dfduu.noalias() += ineqConstraint.dfdu.transpose() * data.linearApproximation.dfdu;
  }
  if (stateConstraint && inputConstraint) {
    cost.dfdux.noalias() += ineqConstraint.dfdu.transpose() * data.linearApproximation.dfdx;
  }
}


void expandDual(const SlackDual& slackDual, InteriorPointMethodData& data) {
  data.dualDirection.array() = - (slackDual.dual.array()*data.slackDirection.array()+data.complementary.array())
                                  / slackDual.slack.array();
}


void expandSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                     const SlackDual& slackDual, const vector_t& dx, const vector_t& du,
                     InteriorPointMethodData& data) {
  data.slackDirection = - data.primalResidual;
  data.slackDirection.noalias() -= ineqConstraint.dfdx * dx;
  data.slackDirection.noalias() -= ineqConstraint.dfdu * du;
  expandDual(slackDual, data);
}


void expandSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                     const SlackDual& slackDual, const vector_t& dx,
                     InteriorPointMethodData& data) {
  data.slackDirection = - data.primalResidual;
  data.slackDirection.noalias() -= ineqConstraint.dfdx * dx;
  expandDual(slackDual, data);
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


scalar_t fractionToBoundaryPrimalStepSize(const SlackDual& slackDual, 
                                          const InteriorPointMethodData& data,
                                          scalar_t marginRate) {
  return fractionToBoundaryStepSize(data.dim, slackDual.slack, data.slackDirection, marginRate);
}


scalar_t fractionToBoundaryDualStepSize(const SlackDual& slackDual, 
                                        const InteriorPointMethodData& data,
                                        scalar_t marginRate) {
  return fractionToBoundaryStepSize(data.dim, slackDual.dual, data.dualDirection, marginRate);
}


scalar_t slackLogBarrier(const SlackDual& slackDual) {
  return - slackDual.barrier * slackDual.slack.array().log().sum();
}

}  // namespace ocs2