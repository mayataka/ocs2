#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <gtest/gtest.h>

#include <ocs2_core/NumericTraits.h>
#include <ocs2_core/reference/ModeSchedule.h>

#include <ocs2_sto/ValidModeSchedule.h>
#include <ocs2_sto/constraint/MinimumDwellTimeConstraint.h>

using namespace ocs2;


TEST(testMinimumDwellTimeConstraint, testMinimumDwellTimeConstraint) {
  auto minimumDwellTimeConstraint = std::unique_ptr<MinimumDwellTimeConstraint>(new MinimumDwellTimeConstraint({{0, 0.1}, {1, 0.1}, {2, 0.2}, }));
  const auto preComp = std::unique_ptr<PreComputation>(new PreComputation);

  const scalar_t initTime  = 0.5;
  const scalar_t finalTime = 5.0; 
  const auto modeSchedule0 = ModeSchedule({}, {0});
  EXPECT_EQ(minimumDwellTimeConstraint->getNumConstraints(initTime, finalTime, modeSchedule0, modeSchedule0), 0);
  const auto linApprox0 = minimumDwellTimeConstraint->getLinearApproximation(initTime, finalTime, modeSchedule0, modeSchedule0, *preComp);
  const auto sizeError0 = checkSize(0, 0, 0, linApprox0, "linApprox0");
  EXPECT_TRUE(sizeError0.empty());

  const auto modeSchedule1 = ModeSchedule({1.0}, {0, 1});
  EXPECT_EQ(minimumDwellTimeConstraint->getNumConstraints(initTime, finalTime, modeSchedule1, modeSchedule1), 2);
  const auto linApprox1 = minimumDwellTimeConstraint->getLinearApproximation(initTime, finalTime, modeSchedule1, modeSchedule1, *preComp);
  const auto sizeError1 = checkSize(2, 1, 0, linApprox1, "linApprox1");
  std::cout << sizeError1 << std::endl;
  EXPECT_TRUE(sizeError1.empty());

  const auto modeSchedule2 = ModeSchedule({1.0, 2.0}, {0, 1, 2});
  EXPECT_EQ(minimumDwellTimeConstraint->getNumConstraints(initTime, finalTime, modeSchedule2, modeSchedule2), 3);
  const auto linApprox2 = minimumDwellTimeConstraint->getLinearApproximation(initTime, finalTime, modeSchedule2, modeSchedule2, *preComp);
  const auto sizeError2 = checkSize(3, 2, 0, linApprox2, "linApprox2");
  std::cout << sizeError2 << std::endl;
  EXPECT_TRUE(sizeError2.empty());
}