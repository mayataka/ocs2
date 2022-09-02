#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <gtest/gtest.h>

#include <ocs2_core/NumericTraits.h>
#include <ocs2_core/reference/ModeSchedule.h>

#include <ocs2_sto/ValidModeSchedule.h>
#include <ocs2_sto/cost/QuadraticStoCost.h>

using namespace ocs2;


TEST(testQuadraticStoCost, testQuadraticStoCost) {
  auto quadraticStoCost = std::unique_ptr<QuadraticStoCost>(new QuadraticStoCost(1.0));

  const scalar_t initTime  = 0.5;
  const scalar_t finalTime = 5.0; 

  const auto referenceModeSchedule = ModeSchedule({1.0, 2.0}, {0, 1, 2});
  const auto stoModeSchedule = referenceModeSchedule;
  const auto preComp = std::unique_ptr<PreComputation>(new PreComputation);

  const auto quadApprox = quadraticStoCost->getQuadraticApproximation(initTime, finalTime, referenceModeSchedule, stoModeSchedule, *preComp);
  EXPECT_DOUBLE_EQ(quadApprox.f, 0.0);

  const auto sizeError = checkSize(referenceModeSchedule.eventTimes.size(), 0, quadApprox, "quadraticStoCostQuadApprox");
  EXPECT_TRUE(sizeError.empty());
}