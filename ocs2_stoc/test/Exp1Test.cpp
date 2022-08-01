#include <gtest/gtest.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>

#include <ocs2_core/initialization/DefaultInitializer.h>
#include <ocs2_oc/test/EXP1.h>

#include <ocs2_stoc/STOC.h>
#include <ocs2_sto/cost/QuadraticStoCost.h>
#include <ocs2_sto/constraint/MinimumDwellTimeConstraint.h>

using namespace ocs2;


class EXP1_StateIneqConstraints final : public StateConstraint {
 public:
  EXP1_StateIneqConstraints(const vector_t& xmin, const vector_t& xmax) : StateConstraint(ConstraintOrder::Linear), xmin_(xmin), xmax_(xmax) {}
  ~EXP1_StateIneqConstraints() override = default;

  EXP1_StateIneqConstraints* clone() const override { return new EXP1_StateIneqConstraints(*this); }

  size_t getNumConstraints(scalar_t time) const override { return 4; }

  vector_t getValue(scalar_t t, const vector_t& x, const PreComputation&) const override {
    vector_t e(4);
    e.head(2) = xmin_ - x;
    e.tail(2) = x - xmax_;
    return e;
  }

  VectorFunctionLinearApproximation getLinearApproximation(scalar_t t, const vector_t& x, const PreComputation& preComp) const override {
    VectorFunctionLinearApproximation e;
    e.f = getValue(t, x, preComp);
    e.dfdx = matrix_t::Zero(4, x.size());
    e.dfdx.topLeftCorner(2, 2) = - matrix_t::Identity(2, 2);
    e.dfdx.bottomRightCorner(2, 2) = matrix_t::Identity(2, 2);
    return e;
  }

 private:
  vector_t xmin_, xmax_;
};


class EXP1_StateInputIneqConstraints final : public StateInputConstraint {
 public:
  EXP1_StateInputIneqConstraints(scalar_t umin, scalar_t umax) : StateInputConstraint(ConstraintOrder::Linear), umin_(umin), umax_(umax) {}
  ~EXP1_StateInputIneqConstraints() override = default;

  EXP1_StateInputIneqConstraints* clone() const override { return new EXP1_StateInputIneqConstraints(*this); }

  size_t getNumConstraints(scalar_t time) const override { return 2; }

  vector_t getValue(scalar_t t, const vector_t& x, const vector_t& u, const PreComputation&) const override {
    vector_t e(2);
    e << (umin_ - u.coeff(0)), (u.coeff(0) - umax_);
    return e;
  }

  VectorFunctionLinearApproximation getLinearApproximation(scalar_t t, const vector_t& x, const vector_t& u,
                                                           const PreComputation& preComp) const override {
    VectorFunctionLinearApproximation e;
    e.f = getValue(t, x, u, preComp);
    e.dfdx = matrix_t::Zero(2, x.size());
    e.dfdu = (matrix_t(2, 1) << -1, 1).finished();
    return e;
  }

 private:
  scalar_t umin_, umax_;
};

TEST(Exp1Test, Unconstrained_FixedSwitchingTimes) {
  static constexpr size_t STATE_DIM = 2;
  static constexpr size_t INPUT_DIM = 1;

  stoc::Settings settings;
  settings.numIteration  = 10;
  settings.useFeedbackPolicy = true;
  settings.dt = 0.01;
  settings.printSolverStatus = true;
  settings.printSolverStatistics = true;
  settings.printLinesearch = true;
  settings.nThreads = 4;

  const scalar_array_t initEventTimes{0.2262, 1.0176};
  const std::vector<size_t> modeSequence{0, 1, 2};
  auto referenceManagerPtr = getExp1ReferenceManager(initEventTimes, modeSequence);
  auto problem = createExp1Problem(referenceManagerPtr);
  auto ipmProblem = ipm::OptimalControlProblem(problem);

  const scalar_t startTime = 0.0;
  const scalar_t finalTime = 3.0;
  const vector_t initState = (vector_t(STATE_DIM) << 2.0, 3.0).finished();

  auto initializerPtr = std::unique_ptr<Initializer>(new DefaultInitializer(INPUT_DIM));

  STOC stoc(settings, ipmProblem, *initializerPtr);
  stoc.setReferenceManager(referenceManagerPtr);
  stoc.run(startTime, initState, finalTime);
  std::cout << stoc.getBenchmarkingInformation() << std::endl;
  std::cout << stoc.getIpmPerformanceIndeces() << std::endl;

  const scalar_t expectedCost = 5.4399;
  EXPECT_NEAR(stoc.getIpmPerformanceIndeces().cost, expectedCost, 0.05); 
  // The error comes from the discretization 

  const auto primalSolution = stoc.primalSolution(finalTime);
  std::cout << "Optimal unconstrained trajectory" << std::endl; 
  for (int i = 0; i < primalSolution.timeTrajectory_.size(); i++) {
    std::cout << "time: " << std::setprecision(4) << primalSolution.timeTrajectory_[i] << "\t state: " << primalSolution.stateTrajectory_[i].transpose()
              << "\t input: " << primalSolution.inputTrajectory_[i].transpose() << std::endl;
  }
}

TEST(Exp1Test, Constrained_FixedSwitchingTimes) {
  static constexpr size_t STATE_DIM = 2;
  static constexpr size_t INPUT_DIM = 1;

  stoc::Settings settings;
  settings.numIteration  = 100;
  settings.primalFeasTol = 1.0e-06;
  settings.dualFeasTol   = 1.0e-06;
  settings.initialBarrierParameter = 1.0e-02;
  settings.targetBarrierParameter = 1.0e-04;
  settings.barrierLinearDecreaseFactor = 0.2;
  settings.barrierSuperlinearDecreasePower = 1.5;
  settings.fractionToBoundaryMargin = 0.995;
  settings.useFeedbackPolicy = true;
  settings.dt = 0.01;
  settings.printSolverStatus = true;
  settings.printSolverStatistics = true;
  settings.printLinesearch = true;
  settings.nThreads = 4;

  const scalar_array_t initEventTimes{0.2262, 1.0176};
  const std::vector<size_t> modeSequence{0, 1, 2};
  auto referenceManagerPtr = getExp1ReferenceManager(initEventTimes, modeSequence);
  auto problem = createExp1Problem(referenceManagerPtr);
  auto ipmProblem = ipm::OptimalControlProblem(problem);

  // add inequality constraints
  const scalar_t umin = -1.0; const scalar_t umax = 1.0;
  std::unique_ptr<StateInputConstraint> stateInputIneqConstraint(new EXP1_StateInputIneqConstraints(umin, umax));
  ipmProblem.inequalityConstraintPtr->add("ubound", std::move(stateInputIneqConstraint));
  const vector_t xmin = (vector_t(2) << -0.0, -0.0).finished(); 
  const vector_t xmax = (vector_t(2) <<  3.0,  4.0).finished(); 
  std::unique_ptr<StateConstraint> stateIneqConstraint(new EXP1_StateIneqConstraints(xmin, xmax));
  std::unique_ptr<StateConstraint> finalStateIneqConstraint(new EXP1_StateIneqConstraints(xmin, xmax));
  ipmProblem.stateInequalityConstraintPtr->add("xbound", std::move(stateIneqConstraint));
  ipmProblem.finalInequalityConstraintPtr->add("xbound", std::move(finalStateIneqConstraint));

  const scalar_t startTime = 0.0;
  const scalar_t finalTime = 3.0;
  const vector_t initState = (vector_t(STATE_DIM) << 2.0, 3.0).finished();

  auto initializerPtr = std::unique_ptr<Initializer>(new DefaultInitializer(INPUT_DIM));

  STOC stoc(settings, ipmProblem, *initializerPtr);
  stoc.setReferenceManager(referenceManagerPtr);
  stoc.run(startTime, initState, finalTime);
  std::cout << stoc.getBenchmarkingInformation() << std::endl;
  std::cout << stoc.getIpmPerformanceIndeces() << std::endl;

  const auto primalSolution = stoc.primalSolution(finalTime);
  std::cout << "Optimal trajectory subject to constraints" << std::endl; 
  std::cout << "xmin: [" << xmin.transpose() << "],  xmax: [" << xmax.transpose() << "],  umin: " << umin << ",  umax: " << umax << std::endl;
  for (int i = 0; i < primalSolution.timeTrajectory_.size(); i++) {
    std::cout << "time: " << std::setprecision(4) << primalSolution.timeTrajectory_[i] << "\t state: " << primalSolution.stateTrajectory_[i].transpose()
              << "\t input: " << primalSolution.inputTrajectory_[i].transpose() << std::endl;
  }
  // check constraint satisfaction
  for (const auto& e : primalSolution.stateTrajectory_) {
    if (e.size() > 0) {
      EXPECT_TRUE(e.coeff(0) >= xmin.coeff(0));
      EXPECT_TRUE(e.coeff(1) >= xmin.coeff(1));
      EXPECT_TRUE(e.coeff(0) <= xmax.coeff(0));
      EXPECT_TRUE(e.coeff(1) <= xmax.coeff(1));
    }
  }
  for (const auto& e : primalSolution.inputTrajectory_) {
    if (e.size() > 0) {
      EXPECT_TRUE(e.coeff(0) >= umin);
      EXPECT_TRUE(e.coeff(0) <= umax);
    }
  }

  // test reset subject to the constraints
  const vector_t randomInitState = vector_t::Random(STATE_DIM) + 0.5 * (xmin + xmax);
  stoc.reset();
  stoc.setReferenceManager(referenceManagerPtr);
  stoc.run(startTime, randomInitState, finalTime);
}


// TEST(Exp1Test, Unconstrained_SwitchingTimesOptimization) {
//   static constexpr size_t STATE_DIM = 2;
//   static constexpr size_t INPUT_DIM = 1;

//   stoc::Settings settings;
//   settings.numIteration  = 10;
//   settings.useFeedbackPolicy = true;
//   settings.dt = 0.01;
//   settings.printSolverStatus = true;
//   settings.printSolverStatistics = true;
//   settings.printLinesearch = true;
//   settings.nThreads = 4;

//   const scalar_array_t initEventTimes{1.0, 2.0};
//   const std::vector<size_t> modeSequence{0, 1, 2};
//   auto referenceManagerPtr = getExp1ReferenceManager(initEventTimes, modeSequence);
//   auto problem = createExp1Problem(referenceManagerPtr);
//   auto ipmProblem = ipm::OptimalControlProblem(problem);
//   auto quadraticStoCost = std::unique_ptr<QuadraticStoCost>(new QuadraticStoCost(0.0));
//   ipmProblem.stoCostPtr->add("quadraticStoCost", std::move(quadraticStoCost));
//   auto minimumDwellTimeConstraint = std::unique_ptr<MinimumDwellTimeConstraint>(new MinimumDwellTimeConstraint({{0, 0.01}, {1, 0.01}, {2, 0.01}}));
//   ipmProblem.stoConstraintPtr->add("minimumDwellTimeConstraint", std::move(minimumDwellTimeConstraint));

//   const scalar_t startTime = 0.0;
//   const scalar_t finalTime = 3.0;
//   const vector_t initState = (vector_t(STATE_DIM) << 2.0, 3.0).finished();

//   auto initializerPtr = std::unique_ptr<Initializer>(new DefaultInitializer(INPUT_DIM));

//   STOC stoc(settings, ipmProblem, *initializerPtr);
//   stoc.setReferenceManager(referenceManagerPtr);
//   stoc.run(startTime, initState, finalTime);
//   std::cout << stoc.getBenchmarkingInformation() << std::endl;
//   std::cout << stoc.getIpmPerformanceIndeces() << std::endl;

//   const scalar_t expectedCost = 5.4399;
//   EXPECT_NEAR(stoc.getIpmPerformanceIndeces().cost, expectedCost, 0.05); 
//   // The error comes from the discretization 

//   const auto primalSolution = stoc.primalSolution(finalTime);
//   std::cout << "Optimal unconstrained trajectory" << std::endl; 
//   for (int i = 0; i < primalSolution.timeTrajectory_.size(); i++) {
//     std::cout << "time: " << std::setprecision(4) << primalSolution.timeTrajectory_[i] << "\t state: " << primalSolution.stateTrajectory_[i].transpose()
//               << "\t input: " << primalSolution.inputTrajectory_[i].transpose() << std::endl;
//   }
// }

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
