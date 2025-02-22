#include <gtest/gtest.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>

#include <ocs2_core/initialization/DefaultInitializer.h>
#include <ocs2_oc/test/EXP0.h>

#include <ocs2_sto/cost/QuadraticStoCost.h>
#include <ocs2_sto/constraint/MinimumDwellTimeConstraint.h>
#include <ocs2_stoc/STOC.h>

using namespace ocs2;

class EXP0_StateIneqConstraints final : public StateConstraint {
 public:
  EXP0_StateIneqConstraints(const vector_t& xmin, const vector_t& xmax) : StateConstraint(ConstraintOrder::Linear), xmin_(xmin), xmax_(xmax) {}
  ~EXP0_StateIneqConstraints() override = default;

  EXP0_StateIneqConstraints* clone() const override { return new EXP0_StateIneqConstraints(*this); }

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


class EXP0_StateInputIneqConstraints final : public StateInputConstraint {
 public:
  EXP0_StateInputIneqConstraints(scalar_t umin, scalar_t umax) : StateInputConstraint(ConstraintOrder::Linear), umin_(umin), umax_(umax) {}
  ~EXP0_StateInputIneqConstraints() override = default;

  EXP0_StateInputIneqConstraints* clone() const override { return new EXP0_StateInputIneqConstraints(*this); }

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

TEST(Exp0Test, Unconstrained_FixedSwitchingTimes) {
  static constexpr size_t STATE_DIM = 2;
  static constexpr size_t INPUT_DIM = 1;

  stoc::Settings settings;
  settings.numIteration = 100;
  settings.useFeedbackPolicy = true;
  settings.dt = 0.01;
  settings.printSolverStatus = true;
  settings.printSolverStatistics = true;
  settings.printLinesearch = true;
  settings.nThreads = 4;

  const scalar_array_t initEventTimes{0.1897};
  const std::vector<size_t> modeSequence{0, 1};
  auto referenceManagerPtr = getExp0ReferenceManager(initEventTimes, modeSequence);
  auto problem = createExp0Problem(referenceManagerPtr);
  auto ipmProblem = ipm::OptimalControlProblem(problem);

  const scalar_t startTime = 0.0;
  const scalar_t finalTime = 2.0;
  const vector_t initState = (vector_t(STATE_DIM) << 0.0, 2.0).finished();

  auto initializerPtr = std::unique_ptr<Initializer>(new DefaultInitializer(INPUT_DIM));

  STOC solver(settings, ipmProblem, *initializerPtr);
  solver.setReferenceManager(referenceManagerPtr);
  solver.run(startTime, initState, finalTime);
  std::cout << solver.getBenchmarkingInformation() << std::endl;
  std::cout << solver.getIpmPerformanceIndeces() << std::endl;

  const scalar_t expectedCost = 9.766;
  EXPECT_NEAR(solver.getIpmPerformanceIndeces().cost, expectedCost, 0.05); 
  // The error comes from the discretization 
  EXPECT_EQ(solver.getNumIterations(), 2); 
  // Since the subsystems are linear subject to no constraints, converges only with a Gauss-Newton iteration for arbitrary initial state.

  const auto primalSolution = solver.primalSolution(finalTime);
  std::cout << "Optimal unconstrained trajectory" << std::endl; 
  for (int i = 0; i < primalSolution.timeTrajectory_.size(); i++) {
    std::cout << "time: " << std::setprecision(4) << primalSolution.timeTrajectory_[i] << "\t state: " << primalSolution.stateTrajectory_[i].transpose()
              << "\t input: " << primalSolution.inputTrajectory_[i].transpose() << std::endl;
  }

  const vector_t randomInitState = vector_t::Random(STATE_DIM);
  solver.reset();
  solver.setReferenceManager(referenceManagerPtr);
  solver.run(startTime, randomInitState, finalTime);
  EXPECT_EQ(solver.getNumIterations(), 2); 
  // Since the subsystems are linear subject to no constraints, converges only with a Gauss-Newton iteration for arbitrary initial state.
}

TEST(Exp0Test, Constrained_FixedSwitchingTimes) {
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

  const scalar_array_t initEventTimes{0.1897};
  const std::vector<size_t> modeSequence{0, 1};
  auto referenceManagerPtr = getExp0ReferenceManager(initEventTimes, modeSequence);
  auto problem = createExp0Problem(referenceManagerPtr);
  auto ipmProblem = ipm::OptimalControlProblem(problem);

  // add inequality constraints
  const scalar_t umin = -7.5; const scalar_t umax = 7.5;
  std::unique_ptr<StateInputConstraint> stateInputIneqConstraint(new EXP0_StateInputIneqConstraints(umin, umax));
  ipmProblem.inequalityConstraintPtr->add("ubound", std::move(stateInputIneqConstraint));
  const vector_t xmin = (vector_t(2) << -7.5, -7.5).finished(); 
  const vector_t xmax = (vector_t(2) <<  7.5,  7.5).finished(); 
  std::unique_ptr<StateConstraint> stateIneqConstraint(new EXP0_StateIneqConstraints(xmin, xmax));
  std::unique_ptr<StateConstraint> finalStateIneqConstraint(new EXP0_StateIneqConstraints(xmin, xmax));
  ipmProblem.stateInequalityConstraintPtr->add("xbound", std::move(stateIneqConstraint));
  ipmProblem.finalInequalityConstraintPtr->add("xbound", std::move(finalStateIneqConstraint));

  const scalar_t startTime = 0.0;
  const scalar_t finalTime = 2.0;
  const vector_t initState = (vector_t(STATE_DIM) << 0.0, 2.0).finished();

  auto initializerPtr = std::unique_ptr<Initializer>(new DefaultInitializer(INPUT_DIM));

  STOC solver(settings, ipmProblem, *initializerPtr);
  solver.setReferenceManager(referenceManagerPtr);
  solver.run(startTime, initState, finalTime);
  std::cout << solver.getBenchmarkingInformation() << std::endl;
  std::cout << solver.getIpmPerformanceIndeces() << std::endl;

  const auto primalSolution = solver.primalSolution(finalTime);
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
  solver.reset();
  solver.setReferenceManager(referenceManagerPtr);
  solver.run(startTime, randomInitState, finalTime);
}

TEST(Exp0Test, Unconstrained_SwitchingTimeOptimization) {
  static constexpr size_t STATE_DIM = 2;
  static constexpr size_t INPUT_DIM = 1;

  stoc::Settings settings;
  settings.numIteration = 100;
  settings.useFeedbackPolicy = true;
  settings.dt = 0.01;
  settings.printSolverStatus = true;
  settings.printSolverStatistics = true;
  settings.printLinesearch = true;
  settings.printSwitchingTimeOptimization = true;
  settings.nThreads = 4;
  settings.initialBarrierParameter = 1.0e-02;
  settings.targetBarrierParameter = 1.0e-04;

  settings.stoEnabledModeSwitches = {{0, 1}, };
  settings.maxTimeInterval = 0.015;
  settings.meshRefinementPrimalFeasTol = 1.0e-02; 
  settings.meshRefinementDualFeasTol = 1.0e-02; 

  settings.switchingTimeTrustRegionRadius = 0.1;
  settings.enableSwitchingTimeTrustRegion = true; 

  settings.printSwitchingTimeOptimization = true; 

  const scalar_array_t initEventTimes{1.0};
  const std::vector<size_t> modeSequence{0, 1};
  auto internalReferenceManagerPtr = getExp0ReferenceManager(initEventTimes, modeSequence);
  auto problem = createExp0Problem(internalReferenceManagerPtr);
  auto ipmProblem = ipm::OptimalControlProblem(problem);
  auto quadraticStoCost = std::unique_ptr<QuadraticStoCost>(new QuadraticStoCost(0.0));
  ipmProblem.stoCostPtr->add("quadraticStoCost", std::move(quadraticStoCost));
  auto minimumDwellTimeConstraint = std::unique_ptr<MinimumDwellTimeConstraint>(new MinimumDwellTimeConstraint({{0, 0.01}, {1, 0.01}}));
  ipmProblem.stoConstraintPtr->add("minimumDwellTimeConstraint", std::move(minimumDwellTimeConstraint));

  const scalar_t startTime = 0.0;
  const scalar_t finalTime = 2.0;
  const vector_t initState = (vector_t(STATE_DIM) << 0.0, 2.0).finished();

  auto initializerPtr = std::unique_ptr<Initializer>(new DefaultInitializer(INPUT_DIM));

  STOC solver(settings, ipmProblem, *initializerPtr);
  auto referenceManagerPtr = getExp0ReferenceManager(initEventTimes, modeSequence);
  solver.setReferenceManager(referenceManagerPtr);
  solver.setInternalReferenceManager(internalReferenceManagerPtr);
  solver.run(startTime, initState, finalTime);
  std::cout << solver.getBenchmarkingInformation() << std::endl;
  std::cout << solver.getIpmPerformanceIndeces() << std::endl;
  std::cout << "\n========= Optimized modeSchedule: =========\n" << solver.getReferenceManager().getModeSchedule() << "\n" << std::endl;

  const scalar_array_t expectedEventTimes = {0.1897};
  EXPECT_NEAR(solver.getReferenceManager().getModeSchedule().eventTimes[0], expectedEventTimes[0], 0.05); 

  const scalar_t expectedCost = 9.766;
  EXPECT_NEAR(solver.getIpmPerformanceIndeces().cost, expectedCost, 0.05); 
}

TEST(Exp0Test, Constrained_SwitchingTimeOptimization) {
  static constexpr size_t STATE_DIM = 2;
  static constexpr size_t INPUT_DIM = 1;

  stoc::Settings settings;
  settings.numIteration = 100;
  settings.useFeedbackPolicy = true;
  settings.dt = 0.01;
  settings.printSolverStatus = true;
  settings.printSolverStatistics = true;
  settings.printLinesearch = true;
  settings.printSwitchingTimeOptimization = true;
  settings.nThreads = 4;
  settings.initialBarrierParameter = 1.0e-02;
  settings.targetBarrierParameter = 1.0e-04;

  settings.stoEnabledModeSwitches = {{0, 1}, };
  settings.maxTimeInterval = 0.015;
  settings.meshRefinementPrimalFeasTol = 1.0e-02; 
  settings.meshRefinementDualFeasTol = 1.0e-02; 

  settings.switchingTimeTrustRegionRadius = 0.1;
  settings.enableSwitchingTimeTrustRegion = true; 

  settings.printSwitchingTimeOptimization = true; 

  const scalar_array_t initEventTimes{0.1897};
  const std::vector<size_t> modeSequence{0, 1};
  auto internalReferenceManagerPtr = getExp0ReferenceManager(initEventTimes, modeSequence);
  auto problem = createExp0Problem(internalReferenceManagerPtr);
  auto ipmProblem = ipm::OptimalControlProblem(problem);
  auto quadraticStoCost = std::unique_ptr<QuadraticStoCost>(new QuadraticStoCost(0.0));
  ipmProblem.stoCostPtr->add("quadraticStoCost", std::move(quadraticStoCost));
  auto minimumDwellTimeConstraint = std::unique_ptr<MinimumDwellTimeConstraint>(new MinimumDwellTimeConstraint({{0, 0.01}, {1, 0.01}}));
  ipmProblem.stoConstraintPtr->add("minimumDwellTimeConstraint", std::move(minimumDwellTimeConstraint));

  // add inequality constraints
  const scalar_t umin = -7.5; const scalar_t umax = 7.5;
  std::unique_ptr<StateInputConstraint> stateInputIneqConstraint(new EXP0_StateInputIneqConstraints(umin, umax));
  ipmProblem.inequalityConstraintPtr->add("ubound", std::move(stateInputIneqConstraint));
  const vector_t xmin = (vector_t(2) << -7.5, -7.5).finished(); 
  const vector_t xmax = (vector_t(2) <<  7.5,  7.5).finished(); 
  std::unique_ptr<StateConstraint> stateIneqConstraint(new EXP0_StateIneqConstraints(xmin, xmax));
  std::unique_ptr<StateConstraint> finalStateIneqConstraint(new EXP0_StateIneqConstraints(xmin, xmax));
  ipmProblem.stateInequalityConstraintPtr->add("xbound", std::move(stateIneqConstraint));
  ipmProblem.finalInequalityConstraintPtr->add("xbound", std::move(finalStateIneqConstraint));

  const scalar_t startTime = 0.0;
  const scalar_t finalTime = 2.0;
  const vector_t initState = (vector_t(STATE_DIM) << 0.0, 2.0).finished();

  auto initializerPtr = std::unique_ptr<Initializer>(new DefaultInitializer(INPUT_DIM));

  STOC solver(settings, ipmProblem, *initializerPtr);
  auto referenceManagerPtr = getExp0ReferenceManager(initEventTimes, modeSequence);
  solver.setReferenceManager(referenceManagerPtr);
  solver.setInternalReferenceManager(internalReferenceManagerPtr);
  solver.run(startTime, initState, finalTime);
  std::cout << solver.getBenchmarkingInformation() << std::endl;
  std::cout << solver.getIpmPerformanceIndeces() << std::endl;
  std::cout << "\n========= Optimized modeSchedule: =========\n" << solver.getReferenceManager().getModeSchedule() << "\n" << std::endl;

  EXPECT_TRUE(solver.getNumIterations() < settings.numIteration);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
