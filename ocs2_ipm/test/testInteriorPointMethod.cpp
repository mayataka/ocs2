#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <gtest/gtest.h>

#include <ocs2_core/Types.h>
#include <ocs2_ipm/ipm/SlackDual.h>
#include <ocs2_ipm/ipm/InteriorPointMethod.h>

using namespace ocs2;



TEST(testInteriorPointMethod, testInteriorPointMethod) {
  const size_t nc = 5;
  const size_t nx = 10;
  const size_t nu = 4;
  const scalar_t barrier = 0.001;
  ipm::SlackDual slackDual;
  slackDual.slack.resize(nc);
  slackDual.dual.resize(nc);
  VectorFunctionLinearApproximation ineqConstraint;
  ineqConstraint.setZero(nc, nx, nu);
  ineqConstraint.f.setRandom();
  ipm::initSlackDual(ineqConstraint, slackDual, barrier);
  EXPECT_TRUE(slackDual.slack.minCoeff() > 0);
  EXPECT_TRUE(slackDual.dual.minCoeff() > 0);
}
