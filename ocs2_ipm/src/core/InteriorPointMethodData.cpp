#include "ocs2_ipm/core/InteriorPointMethodData.h"

#include <iostream>

namespace ocs2 {
namespace ipm {

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
  complementarySlackness.resize(nc);
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
  complementarySlackness.setZero(nc);
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
  if (data.complementarySlackness.size() != constraintDim) {
    errorDescription << dataName << ".complementarySlackness.size() != " << constraintDim << "\n";
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
std::string checkSize(int constraintDim, int stateDim, int inputDim, 
                      const InteriorPointMethodData& data, const std::string& dataName) {
  std::stringstream errorDescription;

  errorDescription << checkSize(constraintDim, data, dataName);
  errorDescription << checkSize(constraintDim, stateDim, inputDim, data.linearApproximation, dataName+".linearApproximation");

  return errorDescription.str();
}

std::ostream& operator<<(std::ostream& out, const InteriorPointMethodData& data) {
  out << "dim: " << data.dim << '\n';
  out << "costBarrier: " << data.costBarrier << '\n';
  out << "primalResidual: " << data.primalResidual.transpose() << '\n';
  out << "complementarySlackness: " << data.complementarySlackness.transpose() << '\n';
  return out;
}

}  // namespace ipm
}  // namespace ocs2