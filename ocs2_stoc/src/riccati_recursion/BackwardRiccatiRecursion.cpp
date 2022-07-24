#include <ocs2_stoc/riccati_recursion/BackwardRiccatiRecursion.h>

#include <cmath>
#include <exception>

namespace ocs2 {
namespace stoc {

BackwardRiccatiRecursion::BackwardRiccatiRecursion(RiccatiSolverMode riccatiSolverMode, scalar_t dts0_max) 
  : riccatiSolverMode_(riccatiSolverMode),
    dts0_max_(),
    sgm_eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    AtP_(),
    BtP_(),
    GK_(),
    Pf_(),
    Ginv_4_(matrix_t::Zero(4, 4)),
    Ginv_3_(matrix_t::Zero(3, 3)),
    Ginv_2_(matrix_t::Zero(2, 2)),
    Ginv_1_(0),
    llt_(),
    ldlt_() {
  setRegularization(dts0_max);
}


void BackwardRiccatiRecursion::setRegularization(scalar_t dts0_max) {
  if (dts0_max < 0) {
    throw std::out_of_range("dts0_max must be non-negative!");
  }
  dts0_max_ = dts0_max;
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
  const size_t nx = modelData.stateDim;
  const size_t nu = modelData.inputDim;;
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
  switch (nu)
  {
  case 4:
    Ginv_4_ = cost.dfduu.inverse();
    lqrPolicy.K.noalias() = - Ginv_4_ * cost.dfdux;
    lqrPolicy.k.noalias() = - Ginv_4_ * cost.dfdu;
    break;
  case 3:
    Ginv_3_ = cost.dfduu.inverse();
    lqrPolicy.K.noalias() = - Ginv_3_ * cost.dfdux;
    lqrPolicy.k.noalias() = - Ginv_3_ * cost.dfdu;
    break;
  case 2:
    Ginv_2_ = cost.dfduu.inverse();
    lqrPolicy.K.noalias() = - Ginv_2_ * cost.dfdux;
    lqrPolicy.k.noalias() = - Ginv_2_ * cost.dfdu;
    break;
  case 1:
    Ginv_1_ = 1.0 / cost.dfduu.coeff(0, 0);
    lqrPolicy.K.noalias() = - Ginv_1_ * cost.dfdux;
    lqrPolicy.k.noalias() = - Ginv_1_ * cost.dfdu;
    break;
  case 0:
    break;
  default:
    if (riccatiSolverMode_ == RiccatiSolverMode::Speed) {
      llt_.compute(cost.dfduu);
      lqrPolicy.K.noalias() = - llt_.solve(cost.dfdux);
      lqrPolicy.k.noalias() = - llt_.solve(cost.dfdu);
    } else {
      ldlt_.compute(cost.dfduu);
      lqrPolicy.K.noalias() = - ldlt_.solve(cost.dfdux);
      lqrPolicy.k.noalias() = - ldlt_.solve(cost.dfdu);
    }
    break;
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
  riccati.s.noalias() -= cost.dfdux.transpose() * lqrPolicy.k;
}


void BackwardRiccatiRecursion::computeIntermediate(const RiccatiRecursionData& riccatiNext, ipm::ModelData& modelData,
                                                   RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy, const bool sto, const bool stoNext) {
  computeIntermediate(riccatiNext, modelData, riccati, lqrPolicy);
  const auto& dynamics = modelData.dynamics;
  const auto& hamiltonian = modelData.hamiltonian;
  if (sto) {
    const size_t nu = modelData.inputDim;
    riccati.psi_x.noalias()  = AtP_ * hamiltonian.dfdt;
    riccati.psi_x.noalias() += hamiltonian.dhdx;
    riccati.psi_x.noalias() += dynamics.dfdx.transpose() * riccatiNext.Psi;
    if (nu > 0) {
      riccati.psi_u.noalias()  = BtP_ * hamiltonian.dfdt;
      riccati.psi_u.noalias() += hamiltonian.dhdu;
      riccati.psi_u.noalias() += dynamics.dfdu.transpose() * riccatiNext.Psi;
    }
    if (stoNext) {
      riccati.phi_x.noalias() = dynamics.dfdx.transpose() * riccatiNext.Phi;
      if (nu > 0) {
        riccati.phi_u.noalias() = dynamics.dfdu.transpose() * riccatiNext.Phi;
      }
    } else {
      riccati.phi_x.setZero();
      if (nu > 0) {
        riccati.phi_u.setZero();
      }
    }
    switch (nu)
    {
    case 3:
      lqrPolicy.T.noalias() = - Ginv_3_ * riccati.psi_u;
      if (stoNext) {
        lqrPolicy.W.noalias() = - Ginv_3_ * riccati.phi_u;
      } else {
        lqrPolicy.W.setZero();
      }
      break;
    case 2:
      lqrPolicy.T.noalias() = - Ginv_2_ * riccati.psi_u;
      if (stoNext) {
        lqrPolicy.W.noalias() = - Ginv_2_ * riccati.phi_u;
      } else {
        lqrPolicy.W.setZero();
      }
      break;
    case 1:
      lqrPolicy.T.noalias() = - Ginv_1_ * riccati.psi_u;
      if (stoNext) {
        lqrPolicy.W.noalias() = - Ginv_1_ * riccati.phi_u;
      } else {
        lqrPolicy.W.setZero();
      }
      break;
    case 0:
      break;
    default:
      if (riccatiSolverMode_ == RiccatiSolverMode::Speed) {
        lqrPolicy.T.noalias() = - llt_.solve(riccati.psi_u);
        if (stoNext) {
          lqrPolicy.W.noalias() = - llt_.solve(riccati.phi_u);
        } else {
          lqrPolicy.W.setZero();
        }
      } else {
        lqrPolicy.T.noalias() = - ldlt_.solve(riccati.psi_u);
        if (stoNext) {
          lqrPolicy.W.noalias() = - ldlt_.solve(riccati.phi_u);
        } else {
          lqrPolicy.W.setZero();
        }
      }
      break;
    }
    // The cost-to-go w.r.t. switching times
    // Vtx
    riccati.Psi              = riccati.psi_x;
    if (nu > 0) {
      riccati.Psi.noalias() += lqrPolicy.K.transpose() * riccati.psi_u;
    }
    if (stoNext) {
      riccati.Phi             = riccati.phi_x;
      if (nu > 0) {
        riccati.Phi.noalias() += lqrPolicy.K.transpose() * riccati.phi_u;
      }
    } else {
      riccati.Phi.setZero();
    }
    // Vtt
    Pf_.noalias() = riccatiNext.P * hamiltonian.dfdt;
    riccati.xi    = hamiltonian.dfdt.dot(Pf_);
    riccati.xi   += hamiltonian.dhdt; 
    riccati.xi   += 2 * riccatiNext.Psi.dot(hamiltonian.dfdt);
    if (nu > 0) {
      riccati.xi   += lqrPolicy.T.dot(riccati.psi_u);
    }
    riccati.xi   += riccatiNext.xi;
    if (stoNext) {
      riccati.chi  = hamiltonian.dhdt;
      // TODO: change this to 
      // riccati.chi  = cost.dfdtt_prev;
      riccati.chi += riccatiNext.Phi.dot(hamiltonian.dfdt);
      riccati.chi += lqrPolicy.T.dot(riccati.phi_u);
      riccati.chi += riccatiNext.chi;
      riccati.rho  = lqrPolicy.W.dot(riccati.phi_u);
      riccati.rho += riccatiNext.rho;
    } else {
      riccati.chi = 0.0;
      riccati.rho = 0.0;
    }
    // Vt
    Pf_.noalias() = riccatiNext.P * dynamics.f - riccatiNext.s;
    riccati.eta   = hamiltonian.dfdt.dot(Pf_);
    riccati.eta  += hamiltonian.h;
    riccati.eta  += riccatiNext.Psi.dot(dynamics.f);
    riccati.eta  += riccati.psi_u.dot(lqrPolicy.k);
    riccati.eta  += riccatiNext.eta;
    if (stoNext) {
      riccati.iota  = riccatiNext.Phi.dot(dynamics.f);
      riccati.iota += riccati.phi_u.dot(lqrPolicy.k);
      riccati.iota += riccatiNext.iota;
    } else {
      riccati.iota = 0.0;
    }
  } else {
    riccati.Psi.setZero();
    riccati.xi = 0.;
    riccati.chi = 0.;
    riccati.eta = 0.;
  }
}

void BackwardRiccatiRecursion::computePreJump(const RiccatiRecursionData& riccatiNext, ipm::ModelData& modelData,
                                              RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy, StoPolicy& stoPolicy, 
                                              const bool sto, const bool stoNext) {
  computeIntermediate(riccatiNext, modelData, riccati, lqrPolicy, sto, stoNext);
  modifyPreJump(riccati, stoPolicy, stoNext);
}

void BackwardRiccatiRecursion::modifyPreJump(RiccatiRecursionData& riccati, StoPolicy& stoPolicy, const bool stoNext) const {
  const size_t nx = riccati.s.size();
  stoPolicy.resize(nx);
  if (stoNext) {
    double sgm = riccati.xi - 2.0 * riccati.chi + riccati.rho;
    if ((sgm*dts0_max_) < std::abs(riccati.eta-riccati.iota) || sgm < sgm_eps_) {
      sgm = std::abs(sgm) + std::abs(riccati.eta-riccati.iota) / dts0_max_;
    }
    stoPolicy.dtsdx  = - (1.0/sgm) * (riccati.Psi-riccati.Phi);
    stoPolicy.dtsdts =   (1.0/sgm) * (riccati.xi-riccati.chi);
    stoPolicy.dts0   = - (1.0/sgm) * (riccati.eta-riccati.iota);
    riccati.s.noalias()   += (1.0/sgm) * (riccati.Psi-riccati.Phi) * (riccati.eta-riccati.iota);
    riccati.Phi.noalias() -= (1.0/sgm) * (riccati.Psi-riccati.Phi) * (riccati.xi-riccati.chi);
    riccati.rho  = riccati.xi - (1.0/sgm) * (riccati.xi-riccati.chi) * (riccati.xi-riccati.chi);
    riccati.iota = riccati.eta - (1.0/sgm) * (riccati.xi-riccati.chi)  * (riccati.eta-riccati.iota);
  }
  riccati.Phi = riccati.Psi;
  riccati.Psi.setZero();
  riccati.rho = riccati.xi;
  riccati.xi = 0.0;
  riccati.chi = 0.0;
  riccati.iota = riccati.eta;
  riccati.eta = 0.0;
}

} // namespace stoc
} // namespace ocs2