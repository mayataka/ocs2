#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <gtest/gtest.h>

#include <ocs2_core/Types.h>

#include <ocs2_ipm/core/SlackDual.h>
#include <ocs2_ipm/core/SlackDualDirection.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>

using namespace ocs2;


TEST(testInteriorPointMethod, testInteriorPointMethod) {
  const size_t nc = 5;
  const size_t nx = 10;
  const size_t nu = 4;
  const scalar_t barrier = 0.001;
  ipm::SlackDual slackDual;
  VectorFunctionLinearApproximation ineqConstraint;
  ineqConstraint.setZero(nc, nx, nu);
  ineqConstraint.f.setRandom();
  ipm::initSlackDual(ineqConstraint, slackDual, barrier);
  EXPECT_TRUE(ipm::checkSize(nc, slackDual, "").empty());
  EXPECT_TRUE(ipm::checkPositive(slackDual, "").empty());
  EXPECT_TRUE((slackDual.slack.array()*slackDual.dual.array()).matrix().isApprox(vector_t::Constant(nc, barrier)));

  ipm::InteriorPointMethodData ipmData;
  ipm::initInteriorPointMethodData(ineqConstraint, ipmData);
  EXPECT_TRUE(ipm::checkSize(nc, ipmData, "").empty());
  EXPECT_TRUE(ipm::checkSize(nc, nx, nu, ipmData, "").empty());

  ScalarFunctionQuadraticApproximation cost;
  cost.setZero(nx, nu);
  ipm::evalPerturbedResidual(ineqConstraint, slackDual, ipmData, cost, true, true);
  EXPECT_TRUE((ineqConstraint.f+slackDual.slack-ipmData.primalResidual).isZero());

  ipm::SlackDualDirection slackDualDirection;
  slackDualDirection.slackDirection.setZero(nc);
  slackDualDirection.dualDirection.setZero(nc);

  EXPECT_NO_THROW(
    ipm::condenseSlackDual(ineqConstraint, slackDual, ipmData, cost, true, true);
    const vector_t dx = vector_t::Random(nx);
    const vector_t du = vector_t::Random(nu);
    ipm::expandSlackDual(ineqConstraint, slackDual, ipmData, dx, du, slackDualDirection);
  );
  const scalar_t primalStepSize = ipm::fractionToBoundaryPrimalStepSize(slackDual, slackDualDirection, ipmData);
  const scalar_t dualStepSize = ipm::fractionToBoundaryDualStepSize(slackDual, slackDualDirection, ipmData);
  EXPECT_TRUE(0.0 <= primalStepSize);
  EXPECT_TRUE(primalStepSize <= 1.0);
  EXPECT_TRUE(0.0 <= dualStepSize);
  EXPECT_TRUE(dualStepSize <= 1.0);
}


TEST(testInteriorPointMethod, testFractionToBoundaryStepSize) {
  const size_t nc = 10;
  const vector_t v = vector_t::Random(nc).cwiseAbs();
  const vector_t dv = vector_t::Random(nc);
  const scalar_t stepSize = ipm::fractionToBoundaryStepSize(nc, v, dv);
  EXPECT_TRUE(0.0 <= stepSize);
  EXPECT_TRUE(stepSize <= 1.0);

  // If the var and dvar are positive, the step size has to be 1.0.
  const vector_t dv1 = vector_t::Random(nc).cwiseAbs();
  const scalar_t stepSize1 = ipm::fractionToBoundaryStepSize(nc, v, dv1);
  EXPECT_DOUBLE_EQ(stepSize1, 1.0);

  // If the var + dvar are negative, the step size has to be less than 1.0.
  const vector_t dv2 = - v - vector_t::Random(nc).cwiseAbs();
  const scalar_t stepSize2 = ipm::fractionToBoundaryStepSize(nc, v, dv2);
  EXPECT_TRUE(0.0 <= stepSize2);
  EXPECT_TRUE(stepSize2 < 1.0);
}