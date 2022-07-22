#include <gtest/gtest.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>

#include <ocs2_core/initialization/DefaultInitializer.h>
#include <ocs2_oc/test/circular_kinematics.h>

#include <ocs2_stoc/STOC.h>

using namespace ocs2;


class CircleKinematics_StateIneqConstraints final : public StateConstraint {
 public:
  CircleKinematics_StateIneqConstraints(const vector_t& xmin, const vector_t& xmax) : StateConstraint(ConstraintOrder::Linear), xmin_(xmin), xmax_(xmax) {}
  ~CircleKinematics_StateIneqConstraints() override = default;

  CircleKinematics_StateIneqConstraints* clone() const override { return new CircleKinematics_StateIneqConstraints(*this); }

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


class CircleKinematics_StateInputIneqConstraints final : public StateInputConstraint {
 public:
  CircleKinematics_StateInputIneqConstraints(const vector_t& umin, const vector_t& umax) : StateInputConstraint(ConstraintOrder::Linear), 
                                                                                           umin_(umin), umax_(umax) {}
  ~CircleKinematics_StateInputIneqConstraints() override = default;

  CircleKinematics_StateInputIneqConstraints* clone() const override { return new CircleKinematics_StateInputIneqConstraints(*this); }

  size_t getNumConstraints(scalar_t time) const override { return 4; }

  vector_t getValue(scalar_t t, const vector_t& x, const vector_t& u, const PreComputation&) const override {
    vector_t e(4);
    e.head(2) = umin_ - u;
    e.tail(2) = u - umax_;
    return e;
  }

  VectorFunctionLinearApproximation getLinearApproximation(scalar_t t, const vector_t& x, const vector_t& u, const PreComputation& preComp) const override {
    VectorFunctionLinearApproximation e;
    e.f = getValue(t, x, u, preComp);
    e.dfdx = matrix_t::Zero(4, x.size());
    e.dfdu = matrix_t::Zero(4, u.size());
    e.dfdu.topLeftCorner(2, 2) = - matrix_t::Identity(2, 2);
    e.dfdu.bottomRightCorner(2, 2) = matrix_t::Identity(2, 2);
    return e;
  }

 private:
  vector_t umin_, umax_;
};


TEST(test_circular_kinematics, solve_projected_EqConstraints) {
  // optimal control problem
  auto problem = createCircularKinematicsProblem("/tmp/sqp_test_generated");
  auto ipmProblem = ipm::OptimalControlProblem(problem);

  // Initializer
  DefaultInitializer zeroInitializer(2);

  // Solver settings
  stoc::Settings settings;
  settings.numIteration  = 10;
  settings.useFeedbackPolicy = true;
  settings.projectStateInputEqualityConstraints = true;
  settings.dt = 0.01;
  settings.printSolverStatus = true;
  settings.printSolverStatistics = true;
  settings.printLinesearch = true;
  settings.nThreads = 4;

  // Additional problem definitions
  const ocs2::scalar_t startTime = 0.0;
  const ocs2::scalar_t finalTime = 1.0;
  const ocs2::vector_t initState = (vector_t(2) << 1.0, 0.0).finished();  // radius 1.0

  // Solve
  STOC stoc(settings, ipmProblem, zeroInitializer);
  stoc.run(startTime, initState, finalTime);
  std::cout << stoc.getBenchmarkingInformation() << std::endl;
  std::cout << stoc.getIpmPerformanceIndeces() << std::endl;

  // Inspect solution
  const auto primalSolution = stoc.primalSolution(finalTime);
  for (int i = 0; i < primalSolution.timeTrajectory_.size(); i++) {
    std::cout << "time: " << std::setprecision(4) << primalSolution.timeTrajectory_[i] << "\t state: " << primalSolution.stateTrajectory_[i].transpose()
              << "\t input: " << primalSolution.inputTrajectory_[i].transpose() << std::endl;
  }

  // check constraint satisfaction
  for (int i = 0; i < primalSolution.timeTrajectory_.size()-1; i++) {
    if (primalSolution.inputTrajectory_[i].size() > 0) {
      EXPECT_NEAR(primalSolution.stateTrajectory_[i].dot(primalSolution.inputTrajectory_[i]), 0.0, 1.0e-06);
    }
  }

  // Check initial condition
  ASSERT_TRUE(primalSolution.stateTrajectory_.front().isApprox(initState));
  ASSERT_DOUBLE_EQ(primalSolution.timeTrajectory_.front(), startTime);
  ASSERT_DOUBLE_EQ(primalSolution.timeTrajectory_.back(), finalTime);

  // Check constraint satisfaction.
  const auto performance = stoc.getPerformanceIndeces();
  ASSERT_LT(performance.dynamicsViolationSSE, 1e-6);
  ASSERT_LT(performance.equalityConstraintsSSE, 1e-6);

  // Check feedback controller
  for (int i = 0; i < primalSolution.timeTrajectory_.size() - 1; i++) {
    const auto t = primalSolution.timeTrajectory_[i];
    const auto& x = primalSolution.stateTrajectory_[i];
    const auto& u = primalSolution.inputTrajectory_[i];
    // Feed forward part
    ASSERT_TRUE(u.isApprox(primalSolution.controllerPtr_->computeInput(t, x)));
  }
}


TEST(test_circular_kinematics, solve_projected_EqConstraints_WithIneqConstraints) {
  // optimal control problem
  auto problem = createCircularKinematicsProblem("/tmp/sqp_test_generated");
  auto ipmProblem = ipm::OptimalControlProblem(problem);
  // add inequality constraints
  const vector_t umin = (vector_t(2) << -1.0e03, -1.0e03).finished(); // no lower bound
  const vector_t umax = (vector_t(2) <<  0.5,  0.5).finished(); 
  std::unique_ptr<StateInputConstraint> stateInputIneqConstraint(new CircleKinematics_StateInputIneqConstraints(umin, umax));
  ipmProblem.inequalityConstraintPtr->add("ubound", std::move(stateInputIneqConstraint));
  const vector_t xmin = (vector_t(2) << -0.5, -0.5).finished(); 
  const vector_t xmax = (vector_t(2) <<  1.0e03,  1.0e03).finished(); // no upper bound
  std::unique_ptr<StateConstraint> stateIneqConstraint(new CircleKinematics_StateIneqConstraints(xmin, xmax));
  std::unique_ptr<StateConstraint> finalStateIneqConstraint(new CircleKinematics_StateIneqConstraints(xmin, xmax));
  ipmProblem.stateInequalityConstraintPtr->add("xbound", std::move(stateIneqConstraint));
  ipmProblem.finalInequalityConstraintPtr->add("xbound", std::move(finalStateIneqConstraint));

  // Initializer
  DefaultInitializer zeroInitializer(2);

  // Solver settings
  stoc::Settings settings;
  settings.numIteration  = 10;
  settings.primalFeasTol = 1.0e-06;
  settings.dualFeasTol   = 1.0e-06;
  settings.initialBarrierParameter = 1.0e-02;
  settings.targetBarrierParameter = 1.0e-04;
  settings.barrierLinearDecreaseFactor = 0.2;
  settings.barrierSuperlinearDecreasePower = 1.5;
  settings.fractionToBoundaryMargin = 0.995;
  settings.useFeedbackPolicy = true;
  settings.projectStateInputEqualityConstraints = true;
  settings.dt = 0.01;
  settings.printSolverStatus = true;
  settings.printSolverStatistics = true;
  settings.printLinesearch = true;
  settings.nThreads = 4;

  // Additional problem definitions
  const ocs2::scalar_t startTime = 0.0;
  const ocs2::scalar_t finalTime = 1.0;
  const ocs2::vector_t initState = (vector_t(2) << 1.0, 0.0).finished();  // radius 1.0

  // Solve
  STOC stoc(settings, ipmProblem, zeroInitializer);
  stoc.run(startTime, initState, finalTime);
  std::cout << stoc.getBenchmarkingInformation() << std::endl;
  std::cout << stoc.getIpmPerformanceIndeces() << std::endl;

  // Inspect solution
  const auto primalSolution = stoc.primalSolution(finalTime);
  for (int i = 0; i < primalSolution.timeTrajectory_.size(); i++) {
    std::cout << "time: " << std::setprecision(4) << primalSolution.timeTrajectory_[i] << "\t state: " << primalSolution.stateTrajectory_[i].transpose()
              << "\t input: " << primalSolution.inputTrajectory_[i].transpose() << std::endl;
  }

  // check constraint satisfaction
  for (int i = 0; i < primalSolution.timeTrajectory_.size()-1; i++) {
    if (primalSolution.inputTrajectory_[i].size() > 0) {
      EXPECT_NEAR(primalSolution.stateTrajectory_[i].dot(primalSolution.inputTrajectory_[i]), 0.0, 1.0e-06);
    }
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
      EXPECT_TRUE(e.coeff(0) >= umin.coeff(0));
      EXPECT_TRUE(e.coeff(1) >= umin.coeff(1));
      EXPECT_TRUE(e.coeff(0) <= umax.coeff(0));
      EXPECT_TRUE(e.coeff(1) <= umax.coeff(1));
    }
  }

  // Check initial condition
  ASSERT_TRUE(primalSolution.stateTrajectory_.front().isApprox(initState));
  ASSERT_DOUBLE_EQ(primalSolution.timeTrajectory_.front(), startTime);
  ASSERT_DOUBLE_EQ(primalSolution.timeTrajectory_.back(), finalTime);

  // Check constraint satisfaction.
  const auto performance = stoc.getPerformanceIndeces();
  ASSERT_LT(performance.dynamicsViolationSSE, 1e-6);
  ASSERT_LT(performance.equalityConstraintsSSE, 1e-6);

  // Check feedback controller
  for (int i = 0; i < primalSolution.timeTrajectory_.size() - 1; i++) {
    const auto t = primalSolution.timeTrajectory_[i];
    const auto& x = primalSolution.stateTrajectory_[i];
    const auto& u = primalSolution.inputTrajectory_[i];
    // Feed forward part
    ASSERT_TRUE(u.isApprox(primalSolution.controllerPtr_->computeInput(t, x)));
  }
}