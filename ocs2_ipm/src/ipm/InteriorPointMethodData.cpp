#include <iostream>

#include <ocs2_ipm/ipm/InteriorPointMethodData.h>

namespace ocs2 {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
InteriorPointMethodData::InteriorPointMethodData(size_t nc, size_t nx, size_t nu) {
  resize(nc, nx, nu);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
InteriorPointMethodData& InteriorPointMethodData::resize(size_t nc, size_t nx, size_t nu) {
  dim = nc;
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

}  // namespace ocs2