#include <iostream>
#include <type_traits>

#include <gtest/gtest.h>

#include <ocs2_ipm/model_data/ModelDataLinearInterpolation.h>

using namespace ocs2;

TEST(testModelData, testModelDataLinearInterpolation) {
  // create data
  const size_t N = 10;
  std::vector<double> timeArray(N);
  std::vector<ipm::ModelData> modelDataBaseArray(N);

  for (size_t i = 0; i < N; i++) {
    double t = 2.0 * i;
    timeArray[i] = t;
    modelDataBaseArray[i].time = t;
    modelDataBaseArray[i].dynamics.f = Eigen::Vector3d::Ones() * t;
    modelDataBaseArray[i].dynamics.dfdx = Eigen::Matrix3d::Ones() * t;
  }

  double time = 5.0;
  // get (index, alpha) pair
  const auto indexAlpha = LinearInterpolation::timeSegment(time, timeArray);

  const scalar_t enquiryScalar = LinearInterpolation::interpolate(indexAlpha, modelDataBaseArray, ipm::model_data::time);
  const vector_t enquiryVector = LinearInterpolation::interpolate(indexAlpha, modelDataBaseArray, ipm::model_data::dynamics_f);
  const matrix_t enquiryMatrix = LinearInterpolation::interpolate(indexAlpha, modelDataBaseArray, ipm::model_data::dynamics_dfdx);

  EXPECT_TRUE(enquiryScalar == time);
  EXPECT_TRUE(enquiryVector.isApprox(Eigen::Vector3d::Ones() * time));
  EXPECT_TRUE(enquiryMatrix.isApprox(Eigen::Matrix3d::Ones() * time));
}

TEST(testModelData, testMovableCopyable) {
  ASSERT_TRUE(std::is_copy_constructible<ipm::ModelData>::value);
  ASSERT_TRUE(std::is_move_constructible<ipm::ModelData>::value);
  ASSERT_TRUE(std::is_nothrow_move_constructible<ipm::ModelData>::value);
}
