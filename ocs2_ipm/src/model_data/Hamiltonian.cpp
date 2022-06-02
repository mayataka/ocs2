#include <iostream>

#include "ocs2_ipm/model_data/Hamiltonian.h"

namespace ocs2 {

Hamiltonian::Hamiltonian(size_t nx, size_t nu) {
  resize(nx, nu);
}


Hamiltonian& Hamiltonian::operator+=(const Hamiltonian& rhs) {
  h += rhs.h;
  dhdt += rhs.dhdt;
  dhdx.noalias() += rhs.dhdx;
  dhdu.noalias() += rhs.dhdu;
  dfdt.noalias() += rhs.dfdt;
  return *this;
}


Hamiltonian& Hamiltonian::operator*=(scalar_t scalar)  {
  h *= scalar;
  dhdt *= scalar;
  dhdx.array() *= scalar;
  dhdu.array() *= scalar;
  dfdt.array() *= scalar;
  return *this;
}


Hamiltonian& Hamiltonian::resize(size_t nx, size_t nu) {
  dhdx.resize(nx);
  dhdu.resize(nu);
  dfdt.resize(nx);
  return *this;
}


Hamiltonian& Hamiltonian::setZero(size_t nx, size_t nu) {
  h = 0.0;
  dhdt = 0.0;
  dhdx.setZero(nx);
  dhdu.setZero(nu);
  dfdt.setZero(nx);
  return *this;
}


Hamiltonian Hamiltonian::Zero(size_t nx, size_t nu) {
  Hamiltonian hamiltonian;
  hamiltonian.setZero(nx, nu);
  return hamiltonian;
}


std::ostream& operator<<(std::ostream& out, const Hamiltonian& hamiltonian) {
  out << "h: " << hamiltonian.h << '\n';
  out << "dhdt: " << hamiltonian.dhdt << '\n';
  out << "dhdx: " << hamiltonian.dhdx.transpose() << '\n';
  out << "dhdu: " << hamiltonian.dhdu.transpose() << '\n';
  out << "dfdt: " << hamiltonian.dfdt.transpose() << '\n';
  return out;
}


std::string checkSize(int stateDim, int inputDim, const Hamiltonian& data, const std::string& dataName) {
  std::stringstream errorDescription;

  if (data.dhdx.size() != stateDim) {
    errorDescription << dataName << ".dhdx.size() != " << stateDim << "\n";
  }
  if (data.dhdu.size() != inputDim) {
    errorDescription << dataName << ".dhdu.size() != " << inputDim << "\n";
  }
  if (data.dfdt.size() != stateDim) {
    errorDescription << dataName << ".dfdt.size() != " << stateDim << "\n";
  }

  return errorDescription.str();
}

}  // namespace ocs2
