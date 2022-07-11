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

  ipm::InteriorPointMethodData ipmData;
  ipmData.setZero(nc, nx, nu);
  EXPECT_TRUE(ipm::checkSize(nc, ipmData, "").empty());

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
}