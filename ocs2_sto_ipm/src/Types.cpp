/******************************************************************************
Copyright (c) 2020, Farbod Farshidian. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <ocs2_sto_ipm/Types.h>

#include <ocs2_core/misc/LinearAlgebra.h>

namespace ocs2 {
namespace sto_ipm {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximationWrapper::ScalarFunctionQuadraticApproximationWrapper(size_t nx, size_t nu) {
  resize(nx, nu);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximationWrapper& ScalarFunctionQuadraticApproximationWrapper::operator+=(const ScalarFunctionQuadraticApproximationWrapper& rhs) {
  base += rhs.base;
  dfdt  += rhs.dfdt;
  dfdtt += rhs.dfdtt;
  dfdxt.noalias() += rhs.dfdxt;
  dfdut.noalias() += rhs.dfdut;
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximationWrapper& ScalarFunctionQuadraticApproximationWrapper::operator*=(scalar_t scalar) {
  base *= scalar;
  dfdt  *= scalar;
  dfdtt *= scalar;
  dfdxt.array() *= scalar;
  dfdut.array() *= scalar;
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximationWrapper& ScalarFunctionQuadraticApproximationWrapper::resize(size_t nx, size_t nu) {
  base.resize(nx, nu);
  dfdxt.resize(nx);
  dfdut.resize(nu);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximationWrapper& ScalarFunctionQuadraticApproximationWrapper::setZero(size_t nx, size_t nu) {
  base.setZero(nx, nu);
  dfdt  = 0.0;
  dfdtt = 0.0;
  dfdxt.setZero(nx);
  dfdut.setZero(nu);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximationWrapper ScalarFunctionQuadraticApproximationWrapper::Zero(size_t nx, size_t nu) {
  ScalarFunctionQuadraticApproximationWrapper f;
  f.setZero(nx, nu);
  return f;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkBeingPSD(const ScalarFunctionQuadraticApproximationWrapper& data, const std::string& dataName) {
  return checkBeingPSD(data.base, dataName);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkSize(int stateDim, int inputDim, const ScalarFunctionQuadraticApproximationWrapper& data, const std::string& dataName) {
  std::stringstream errorDescription;

  if (data.dfdxt.size() != stateDim) {
    errorDescription << dataName << ".dfdxt.size() != " << stateDim << "\n";
  }
  if (data.dfdut.size() != inputDim) {
    errorDescription << dataName << ".dfdut.size() != " << inputDim << "\n";
  }

  return checkSize(stateDim, inputDim, data.base, dataName) + errorDescription.str();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::ostream& operator<<(std::ostream& out, const ScalarFunctionQuadraticApproximationWrapper& f) {
  out << f.base;
  out << "dfdt: "  << f.dfdt  << '\n';
  out << "dfdtt: " << f.dfdtt << '\n';
  out << "dfdxt: " << f.dfdxt.transpose() << '\n';
  out << "dfdut: " << f.dfdut.transpose() << '\n';
  return out;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximationWrapper::VectorFunctionLinearApproximationWrapper(size_t nv, size_t nx, size_t nu) {
  resize(nv, nx, nu);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximationWrapper& VectorFunctionLinearApproximationWrapper::resize(size_t nv, size_t nx, size_t nu) {
  base.resize(nv, nx, nu);
  dfdt.resize(nv);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximationWrapper& VectorFunctionLinearApproximationWrapper::setZero(size_t nv, size_t nx, size_t nu) {
  base.setZero(nv, nx, nu);
  dfdt.setZero(nv);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionLinearApproximationWrapper VectorFunctionLinearApproximationWrapper::Zero(size_t nv, size_t nx, size_t nu) {
  VectorFunctionLinearApproximationWrapper f;
  f.setZero(nv, nx, nu);
  return f;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::ostream& operator<<(std::ostream& out, const VectorFunctionLinearApproximationWrapper& f) {
  out << f.base;
  out << "dfdt:\n" << f.dfdt.transpose() << '\n';
  return out;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::string checkSize(int vectorDim, int stateDim, int inputDim, const VectorFunctionLinearApproximationWrapper& data,
                      const std::string& dataName) {
  std::stringstream errorDescription;

  if (data.dfdt.size() != vectorDim) {
    errorDescription << dataName << ".dfdt.size() != " << vectorDim << "\n";
  }

  return checkSize(vectorDim, stateDim, inputDim, data.base, dataName) + errorDescription.str();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionQuadraticApproximationWrapper::VectorFunctionQuadraticApproximationWrapper(size_t nv, size_t nx, size_t nu) {
  resize(nv, nx, nu);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionQuadraticApproximationWrapper& VectorFunctionQuadraticApproximationWrapper::resize(size_t nv, size_t nx, size_t nu) {
  base.resize(nv, nx, nu);
  dfdt.resize(nv);
  dfdtt.resize(nv);
  dfdxt.resize(nv, nx);
  dfdut.resize(nv, nu);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionQuadraticApproximationWrapper& VectorFunctionQuadraticApproximationWrapper::setZero(size_t nv, size_t nx, size_t nu) {
  base.setZero(nv, nx, nu);
  dfdt.setZero(nv);
  dfdtt.setZero(nv);
  dfdxt.setZero(nv, nx);
  dfdut.setZero(nv, nu);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
VectorFunctionQuadraticApproximationWrapper VectorFunctionQuadraticApproximationWrapper::Zero(size_t nv, size_t nx, size_t nu) {
  VectorFunctionQuadraticApproximationWrapper f;
  f.setZero(nv, nx, nu);
  return f;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::ostream& operator<<(std::ostream& out, const VectorFunctionQuadraticApproximationWrapper& f) {
  out << f.base;
  out << "dfdt:\n" << f.dfdt.transpose() << '\n';
  out << "dfdtt:\n" << f.dfdtt.transpose() << '\n';
  out << "dfdxt:\n" << f.dfdxt << '\n';
  out << "dfdut:\n" << f.dfdut << '\n';
  return out;
}

}  // namespace sto_ipm
}  // namespace ocs2
