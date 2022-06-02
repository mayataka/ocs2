#pragma once

#include <utility>
#include <ocs2_core/Types.h>

namespace ocs2 {
namespace stoc {

struct RiccatiRecursionData {
  matrix_t P;
  vector_t s;
  vector_t psi_x;
  vector_t psi_u;
  vector_t Psi;
  vector_t phi_x;
  vector_t phi_u;
  vector_t Phi;
  scalar_t xi;
  scalar_t chi;
  scalar_t rho;
  scalar_t eta;
  scalar_t iota;

  ScalarFunctionQuadraticApproximation toValueFunctionQuadraticApproximation() const {
    ScalarFunctionQuadraticApproximation valueFunction;
    valueFunction.setZero(s.size(), 0);
    valueFunction.dfdxx = P;
    valueFunction.dfdx = - s;
    return valueFunction;
  }

  void resize(const size_t nx, const size_t nu) {
    P.resize(nx, nx);
    s.resize(nx);
    psi_x.resize(nx);
    psi_u.resize(nu);
    Psi.resize(nx);
    phi_x.resize(nx);
    phi_u.resize(nu);
    Phi.resize(nx);
  }

  void setZero() {
    P.setZero();
    s.setZero();
    psi_x.setZero();
    psi_u.setZero();
    Psi.setZero();
    phi_x.setZero();
    phi_u.setZero();
    Phi.setZero();
    xi = 0.0;
    chi = 0.0;
    rho = 0.0;
    eta = 0.0;
    iota = 0.0;
  }

  void swap(RiccatiRecursionData& other) {
    P.swap(other.P);
    s.swap(other.s);
    psi_x.swap(other.psi_x);
    psi_u.swap(other.psi_u);
    Psi.swap(other.Psi);
    phi_x.swap(other.phi_x);
    phi_u.swap(other.phi_u);
    Phi.swap(other.Phi);
    std::swap(xi, other.xi);
    std::swap(chi, other.chi);
    std::swap(rho, other.rho);
    std::swap(eta, other.eta);
    std::swap(iota, other.iota);
  }
};

} // namespace stoc
} // namespace ocs2