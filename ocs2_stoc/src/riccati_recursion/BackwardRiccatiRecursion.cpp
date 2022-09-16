#include "ocs2_stoc/riccati_recursion/BackwardRiccatiRecursion.h"

#include <cmath>
#include <exception>
#include <iostream>

namespace ocs2 {
namespace stoc {

BackwardRiccatiRecursion::BackwardRiccatiRecursion(RiccatiSolverMode riccatiSolverMode, scalar_t switchingTimeTrustRegionRadius,
                                                   bool enableSwitchingTimeTrustRegion) 
  : riccatiSolverMode_(riccatiSolverMode),
    switchingTimeTrustRegionRadius_(switchingTimeTrustRegionRadius),
    enableSwitchingTimeTrustRegion_(enableSwitchingTimeTrustRegion),
    minQuadraticCoeff_(std::sqrt(std::numeric_limits<double>::epsilon())),
    AtP_(),
    BtP_(),
    GK_(),
    Pf_(),
    Ginv_(),
    llt_(),
    ldlt_() {
  setRegularization(switchingTimeTrustRegionRadius, enableSwitchingTimeTrustRegion);
}


void BackwardRiccatiRecursion::setRegularization(scalar_t switchingTimeTrustRegionRadius, bool enableSwitchingTimeTrustRegion) {
  if (switchingTimeTrustRegionRadius < 0) {
    throw std::out_of_range("switchingTimeTrustRegionRadius must be non-negative!");
  }
  switchingTimeTrustRegionRadius_ = switchingTimeTrustRegionRadius;
  enableSwitchingTimeTrustRegion_ = enableSwitchingTimeTrustRegion;
}


void BackwardRiccatiRecursion::resize(size_t nx, size_t nu) {
  AtP_.resize(nx, nx);
  BtP_.resize(nu, nx);
  GK_.resize(nu, nx);
  Pf_.resize(nx); 
}


void BackwardRiccatiRecursion::computeFinal(const ipm::ModelData& modelData, RiccatiRecursionData& riccati) {
  riccati.resize(modelData.stateDim, 0);
  const auto& cost = modelData.cost;
  riccati.P = cost.dfdxx;
  riccati.s = - cost.dfdx;
}


void BackwardRiccatiRecursion::computeIntermediate(const RiccatiRecursionData& riccatiNext, ipm::ModelData& modelData, 
                                                   RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy) {
  auto& cost = modelData.cost;
  const auto& dynamics = modelData.dynamics;
  const auto nx = modelData.stateDim;
  const auto nu = modelData.inputDim;
  riccati.resize(nx, nu);
  lqrPolicy.resize(nx, nu);
  this->resize(nx, nu);
  AtP_.noalias() = dynamics.dfdx.transpose() * riccatiNext.P;
  if (nu > 0) {
    BtP_.noalias() = dynamics.dfdu.transpose() * riccatiNext.P;
  }
  // Factorize F
  cost.dfdxx.noalias() += AtP_ * dynamics.dfdx;
  if (nu > 0) {
    // Factorize H^T
    cost.dfdux.noalias() += BtP_ * dynamics.dfdx;
    // Factorize G
    cost.dfduu.noalias() += BtP_ * dynamics.dfdu;
    // Factorize vector term
    cost.dfdu.noalias() += BtP_ * dynamics.f;
    cost.dfdu.noalias() -= dynamics.dfdu.transpose() * riccatiNext.s; 
  }
  if (nu > 4) {
    if (riccatiSolverMode_ == RiccatiSolverMode::Speed) {
      llt_.compute(cost.dfduu);
      lqrPolicy.K.noalias() = - llt_.solve(cost.dfdux);
      lqrPolicy.k.noalias() = - llt_.solve(cost.dfdu);
    } else {
      ldlt_.compute(cost.dfduu);
      lqrPolicy.K.noalias() = - ldlt_.solve(cost.dfdux);
      lqrPolicy.k.noalias() = - ldlt_.solve(cost.dfdu);
    }
  }
  else if (nu > 0) {
    Ginv_ = cost.dfduu.inverse();
    lqrPolicy.K.noalias() = - Ginv_ * cost.dfdux;
    lqrPolicy.k.noalias() = - Ginv_ * cost.dfdu;
  }
  if (nu > 0) {
    GK_.noalias() = cost.dfduu * lqrPolicy.K; 
    cost.dfdxx.noalias() -= lqrPolicy.K.transpose() * GK_;
  }
  // Riccati factorization matrix with preserving the symmetry
  riccati.P = 0.5 * (cost.dfdxx + cost.dfdxx.transpose());
  // Riccati factorization vector
  riccati.s.noalias()  = dynamics.dfdx.transpose() * riccatiNext.s;
  riccati.s.noalias() -= AtP_ * dynamics.f;
  riccati.s.noalias() -= cost.dfdx;
  if (nu > 0) {
    riccati.s.noalias() -= cost.dfdux.transpose() * lqrPolicy.k;
  }
}


void BackwardRiccatiRecursion::computeIntermediate(const RiccatiRecursionData& riccatiNext, ipm::ModelData& modelData,
                                                   RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy, 
                                                   const bool sto, const bool stoNext) {
  computeIntermediate(riccatiNext, modelData, riccati, lqrPolicy);
  const auto& dynamics = modelData.dynamics;
  const auto& hamiltonian = modelData.hamiltonian;

  const auto nx = modelData.stateDim;
  const auto nu = modelData.inputDim;
  if (!sto) {
    riccati.psi_x.setZero(nx);
    riccati.psi_u.setZero(nu);
    lqrPolicy.T.setZero(nu);
  } else {
    riccati.psi_x.noalias()  = AtP_ * hamiltonian.dfdt;
    riccati.psi_x.noalias() += hamiltonian.dhdx;
    riccati.psi_x.noalias() += dynamics.dfdx.transpose() * riccatiNext.Psi;
    if (nu > 0) {
      riccati.psi_u.noalias()  = BtP_ * hamiltonian.dfdt;
      riccati.psi_u.noalias() += hamiltonian.dhdu;
      riccati.psi_u.noalias() += dynamics.dfdu.transpose() * riccatiNext.Psi;
    }
    if (nu > 4) {
      if (riccatiSolverMode_ == RiccatiSolverMode::Speed) {
        lqrPolicy.T.noalias() = - llt_.solve(riccati.psi_u);
      } else {
        lqrPolicy.T.noalias() = - ldlt_.solve(riccati.psi_u);
      }
    }
    else if (nu > 0) {
      lqrPolicy.T.noalias() = - Ginv_ * riccati.psi_u;
    }
  }
  if (!stoNext) {
    riccati.phi_x.setZero(nx);
    riccati.phi_u.setZero(nu);
    lqrPolicy.W.setZero(nu);
  } else {
    riccati.phi_x.noalias() = dynamics.dfdx.transpose() * riccatiNext.Phi;
    if (nu > 0) {
      riccati.phi_u.noalias() = dynamics.dfdu.transpose() * riccatiNext.Phi;
    }
    if (nu > 4) {
      if (riccatiSolverMode_ == RiccatiSolverMode::Speed) {
        lqrPolicy.W.noalias() = - llt_.solve(riccati.phi_u);
      } else {
        lqrPolicy.W.noalias() = - ldlt_.solve(riccati.phi_u);
      }
    }
    else if (nu > 0) {
      lqrPolicy.W.noalias() = - Ginv_ * riccati.phi_u;
    }
  }
  // The cost-to-go w.r.t. switching times
  // Vtx
  riccati.Psi              = riccati.psi_x;
  if (nu > 0) {
    riccati.Psi.noalias() += lqrPolicy.K.transpose() * riccati.psi_u;
  }
  riccati.Phi              = riccati.phi_x;
  if (nu > 0) {
    riccati.Phi.noalias() += lqrPolicy.K.transpose() * riccati.phi_u;
  }
  // Vtt
  riccati.xi  = 0.;
  if (sto) {
    Pf_.noalias() = riccatiNext.P * hamiltonian.dfdt;
    riccati.xi    = hamiltonian.dfdt.dot(Pf_);
    riccati.xi   += hamiltonian.dhdt; 
    riccati.xi   += 2 * riccatiNext.Psi.dot(hamiltonian.dfdt);
    if (nu > 0) {
      riccati.xi += lqrPolicy.T.dot(riccati.psi_u);
    }
    riccati.xi   += riccatiNext.xi;
  }
  riccati.chi = 0.0;
  riccati.rho = 0.0;
  if (stoNext) {
    riccati.chi    = hamiltonian.dhdt;
    // TODO: change this to 
    // riccati.chi  = cost.dfdtt_prev;
    riccati.chi   += riccatiNext.Phi.dot(hamiltonian.dfdt);
    if (nu > 0) {
      riccati.chi += lqrPolicy.T.dot(riccati.phi_u);
    }
    riccati.chi   += riccatiNext.chi;
    riccati.rho    = riccatiNext.rho;
    if (nu > 0) {
      riccati.rho += lqrPolicy.W.dot(riccati.phi_u);
    }
  }
  // Vt
  Pf_.noalias() = riccatiNext.P * dynamics.f - riccatiNext.s;
  riccati.eta = 0.0;
  if (sto) {
    riccati.eta  = hamiltonian.dfdt.dot(Pf_);
    riccati.eta  += hamiltonian.h;
    riccati.eta  += riccatiNext.Psi.dot(dynamics.f);
    if (nu > 0) {
      riccati.eta  += riccati.psi_u.dot(lqrPolicy.k);
    }
    riccati.eta  += riccatiNext.eta;
  }
  riccati.iota = 0.0;
  if (stoNext) {
    riccati.iota  = riccatiNext.Phi.dot(dynamics.f);
    if (nu > 0) {
      riccati.iota += riccati.phi_u.dot(lqrPolicy.k);
    }
    riccati.iota += riccatiNext.iota;
  }
}


void BackwardRiccatiRecursion::computePreJump(const RiccatiRecursionData& riccatiNext, ipm::ModelData& modelData,
                                              RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy, const bool sto, const bool stoNext) {
  computeIntermediate(riccatiNext, modelData, riccati, lqrPolicy, sto, stoNext);
}


void BackwardRiccatiRecursion::computePostJump(const RiccatiRecursionData& riccatiNext, ipm::ModelData& modelData,
                                               RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy, StoPolicy& stoPolicy, 
                                               const bool sto, const bool stoNext) {
  computeIntermediate(riccatiNext, modelData, riccati, lqrPolicy, sto, stoNext);
  modifyPostJump(riccati, stoPolicy, stoNext);
}


void BackwardRiccatiRecursion::modifyPostJump(RiccatiRecursionData& riccati, StoPolicy& stoPolicy, bool computeStoPolicy) const {
  const size_t nx = riccati.s.size();
  stoPolicy.setZero(nx);
  if (computeStoPolicy) {
    scalar_t sgm = riccati.xi - 2.0 * riccati.chi + riccati.rho;
    if (enableSwitchingTimeTrustRegion_ 
        && ((sgm*switchingTimeTrustRegionRadius_) < std::abs(riccati.eta-riccati.iota) || sgm < minQuadraticCoeff_)) {
      sgm = std::abs(sgm) + std::abs(riccati.eta-riccati.iota) / switchingTimeTrustRegionRadius_;
    }
    stoPolicy.dtsdx  = - (1.0/sgm) * (riccati.Psi-riccati.Phi);
    stoPolicy.dtsdts =   (1.0/sgm) * (riccati.xi-riccati.chi);
    stoPolicy.dts0   = - (1.0/sgm) * (riccati.eta-riccati.iota);
    riccati.s.noalias() += (1.0/sgm) * (riccati.Psi-riccati.Phi) * (riccati.eta-riccati.iota);
    riccati.Phi  = riccati.Psi - (1.0/sgm) * (riccati.Psi-riccati.Phi) * (riccati.xi-riccati.chi);
    riccati.rho  = riccati.xi - (1.0/sgm) * (riccati.xi-riccati.chi) * (riccati.xi-riccati.chi);
    riccati.iota = riccati.eta - (1.0/sgm) * (riccati.xi-riccati.chi)  * (riccati.eta-riccati.iota);
    riccati.Psi.setZero();
    riccati.xi  = 0.0;
    riccati.chi = 0.0;
    riccati.eta = 0.0;
  }
  else {
    riccati.Phi  = riccati.Psi;
    riccati.rho  = riccati.xi;
    riccati.iota = riccati.eta;
    riccati.Psi.setZero();
    riccati.xi   = 0.0;
    riccati.chi  = 0.0;
    riccati.eta  = 0.0;
  }
}

} // namespace stoc
} // namespace ocs2