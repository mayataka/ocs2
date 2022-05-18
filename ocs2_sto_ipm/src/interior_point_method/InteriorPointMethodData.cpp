#include <iostream>

#include <ocs2_sto_ipm/interior_point_method/InteriorPointMethodData.h>

namespace ocs2 {
namespace sto_ipm {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
InteriorPointMethodData::InteriorPointMethodData(size_t nc, size_t nx, size_t nu) {
  resize(nc, nx, nu);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
InteriorPointMethodData::InteriorPointMethodData(size_t nc, size_t nx, size_t nu, scalar_t barrier) {
  setBarrier(nc, nx, nu, barrier);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
InteriorPointMethodData& InteriorPointMethodData::resize(size_t nc, size_t nx, size_t nu) {
  dim = nc;
  slack.resize(nc);
  dual.resize(nc);
  primalResidual.resize(nc);
  complementary.resize(nc);
  slackDirection.resize(nc);
  dualDirection.resize(nc);
  cond.resize(nc);
  dualDivSlack.resize(nc);
  linearApproximation.resize(nc, nx, nu);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
InteriorPointMethodData& InteriorPointMethodData::setZero(size_t nc, size_t nx, size_t nu) {
  dim = nc;
  slack.setZero(nc);
  dual.setZero(nc);
  primalResidual.setZero(nc);
  complementary.setZero(nc);
  slackDirection.setZero(nc);
  dualDirection.setZero(nc);
  cond.setZero(nc);
  dualDivSlack.setZero(nc);
  linearApproximation.setZero(nc, nx, nu);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
InteriorPointMethodData& InteriorPointMethodData::setBarrier(scalar_t inputBarrier) {
  barrier = inputBarrier;
  slack.fill(std::sqrt(barrier));
  dual.fill(std::sqrt(barrier));
  primalResidual.setZero();
  complementary.setZero();
  slackDirection.setZero();
  dualDirection.setZero();
  cond.setZero();
  dualDivSlack.setZero();
  linearApproximation.f.setZero();
  linearApproximation.dfdx.setZero();
  linearApproximation.dfdu.setZero();
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
InteriorPointMethodData& InteriorPointMethodData::setBarrier(size_t nc, size_t nx, size_t nu, scalar_t inputBarrier) {
  resize(nc, nx, nu);
  setBarrier(inputBarrier);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
InteriorPointMethodData InteriorPointMethodData::Zero(size_t nc, size_t nx, size_t nu) {
  InteriorPointMethodData data(nc, nx, nu);
  return data;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkSize(int constraintDim, const InteriorPointMethodData& data, const std::string& dataName) {
  std::stringstream errorDescription;

  if (data.dim != constraintDim) {
    errorDescription << dataName << ".dim != " << constraintDim << "\n";
  }
  if (data.slack.size() != constraintDim) {
    errorDescription << dataName << ".slack.size() != " << constraintDim << "\n";
  }
  if (data.dual.size() != constraintDim) {
    errorDescription << dataName << ".dual.size() != " << constraintDim << "\n";
  }
  if (data.primalResidual.size() != constraintDim) {
    errorDescription << dataName << ".primalResidual.size() != " << constraintDim << "\n";
  }
  if (data.complementary.size() != constraintDim) {
    errorDescription << dataName << ".complementary.size() != " << constraintDim << "\n";
  }
  if (data.slackDirection.size() != constraintDim) {
    errorDescription << dataName << ".slackDirection.size() != " << constraintDim << "\n";
  }
  if (data.dualDirection.size() != constraintDim) {
    errorDescription << dataName << ".dualDirection.size() != " << constraintDim << "\n";
  }
  if (data.cond.size() != constraintDim) {
    errorDescription << dataName << ".cond.size() != " << constraintDim << "\n";
  }
  if (data.dualDivSlack.size() != constraintDim) {
    errorDescription << dataName << ".dualDivSlack.size() != " << constraintDim << "\n";
  }
  return errorDescription.str();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkPositive(const InteriorPointMethodData& data, const std::string& dataName) {
  std::stringstream errorDescription;

  if (data.barrier <= 0.0) {
    errorDescription << dataName << ".barrier = " << data.barrier << "\n";
  }
  if (data.slack.minCoeff() <= 0.0) {
    errorDescription << dataName << ".slack = " << data.slack.transpose() << "\n";
  }
  if (data.dual.minCoeff() <= 0.0) {
    errorDescription << dataName << ".dual = " << data.dual.transpose() << "\n";
  }
  return errorDescription.str();
}

}  // namespace sto_ipm
}  // namespace ocs2