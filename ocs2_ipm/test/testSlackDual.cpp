#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <gtest/gtest.h>

#include <ocs2_core/Types.h>
#include <ocs2_ipm/ipm/SlackDual.h>
#include <ocs2_ipm/ipm/InteriorPointMethod.h>

using namespace ocs2;


TEST(testInteriorPointMethod, testSlackDual) {
  const size_t nc = 5;
  const size_t nx = 10;
  const size_t nu = 4;
  const scalar_t barrier = 0.001;
  ipm::SlackDual slackDual;
  VectorFunctionLinearApproximation ineqConstraint;
  ineqConstraint.setZero(nc, nx, nu);
  ineqConstraint.f.setRandom();
  ipm::initSlackDual(ineqConstraint, slackDual, barrier);
  EXPECT_EQ(slackDual.slack.size(), nc);
  EXPECT_EQ(slackDual.dual.size(), nc);
  EXPECT_TRUE(slackDual.slack.minCoeff() > 0);
  EXPECT_TRUE(slackDual.dual.minCoeff() > 0);
  EXPECT_TRUE(ipm::checkSize(nc, slackDual, "").empty());
  EXPECT_TRUE(ipm::checkPositive(slackDual, "").empty());
}
