#include <ocs2_stoc/riccati_recursion/BackwardRiccatiRecursion.h>

#include <cmath>
#include <exception>

namespace ocs2 {
namespace stoc {

BackwardRiccatiRecursion::BackwardRiccatiRecursion(scalar_t dts0_max) 
  : dts0_max_(),
    sgm_eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    AtP_(),
    BtP_(),
    GK_(),
    Pf_(),
    Ginv_3_(),
    Ginv_2_(),
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


void BackwardRiccatiRecursion::compute(const ScalarFunctionQuadraticApproximation& cost, RiccatiRecursionData& riccati) {
  const size_t nx = cost.dfdx.size();
  const size_t nu = cost.dfdu.size();
  riccati.resize(nx, nu);
  riccati.P = cost.dfdxx;
  riccati.s = - cost.dfdx;
}


void BackwardRiccatiRecursion::compute(const RiccatiRecursionData& riccatiNext, const VectorFunctionLinearApproximation& dynamics, 
                                       ScalarFunctionQuadraticApproximation& cost, RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy) {
  const size_t nx = cost.dfdx.size();
  const size_t nu = cost.dfdu.size();
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
    ldlt_.compute(cost.dfduu);
    lqrPolicy.K.noalias() = - ldlt_.solve(cost.dfdux);
    lqrPolicy.k.noalias() = - ldlt_.solve(cost.dfdu);
    break;
  }
}


void BackwardRiccatiRecursion::compute(const RiccatiRecursionData& riccatiNext, const VectorFunctionLinearApproximation& dynamics, 
                                       ScalarFunctionQuadraticApproximation& cost, ipm::Hamiltonian& hamiltonian, 
                                       RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy, const bool sto, const bool stoNext) {
  compute(riccatiNext, dynamics, cost, riccati, lqrPolicy);
  if (sto) {
    const size_t nu = cost.dfdu.size();
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
    }
    else {
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
      }
      else {
        lqrPolicy.W.setZero();
      }
      break;
    case 2:
      lqrPolicy.T.noalias() = - Ginv_2_ * riccati.psi_u;
      if (stoNext) {
        lqrPolicy.W.noalias() = - Ginv_2_ * riccati.phi_u;
      }
      else {
        lqrPolicy.W.setZero();
      }
      break;
    case 1:
      lqrPolicy.T.noalias() = - Ginv_1_ * riccati.psi_u;
      if (stoNext) {
        lqrPolicy.W.noalias() = - Ginv_1_ * riccati.phi_u;
      }
      else {
        lqrPolicy.W.setZero();
      }
      break;
    case 0:
      break;
    default:
      lqrPolicy.T.noalias() = - ldlt_.solve(riccati.psi_u);
      if (stoNext) {
        lqrPolicy.W.noalias() = - ldlt_.solve(riccati.phi_u);
      }
      else {
        lqrPolicy.W.setZero();
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
    }
    else {
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
    }
    else {
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
    }
    else {
      riccati.iota = 0.0;
    }
  }
  else {
    riccati.Psi.setZero();
    riccati.xi = 0.;
    riccati.chi = 0.;
    riccati.eta = 0.;
  }
}


void BackwardRiccatiRecursion::compute(RiccatiRecursionData& riccati, StoPolicy& stoPolicy, const bool stoNext) const {
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


// void BackwardRiccatiRecursion::postEventToPreEvent(const RiccatiRecursionData& riccatiPostEvent, 
//                                                    RiccatiRecursionData& riccatiPreEvent, 
//                                                    StoPolicy& stoPolicy, const bool stoNext) const {
//   riccatiPreEvent.P = riccatiPostEvent.P;
//   riccatiPreEvent.s = riccatiPostEvent.s;
//   riccatiPreEvent.Psi.setZero();
//   riccatiPreEvent.Phi = riccatiPostEvent.Psi;
//   riccatiPreEvent.xi = 0.0;
//   riccatiPreEvent.chi = 0.0;
//   riccatiPreEvent.rho = riccatiPostEvent.xi;
//   riccatiPreEvent.eta = 0.0;
//   riccatiPreEvent.iota = riccatiPostEvent.eta;
//   if (stoNext) {
//     double sgm = riccatiPostEvent.xi - 2.0 * riccatiPostEvent.chi + riccatiPostEvent.rho;
//     if ((sgm*dts0_max_) < std::abs(riccatiPostEvent.eta-riccatiPostEvent.iota) || sgm < sgm_eps_) {
//       sgm = std::abs(sgm) + std::abs(riccatiPostEvent.eta-riccatiPostEvent.iota) / dts0_max_;
//     }
//     stoPolicy.dtsdx  = - (1.0/sgm) * (riccatiPostEvent.Psi-riccatiPostEvent.Phi);
//     stoPolicy.dtsdts =   (1.0/sgm) * (riccatiPostEvent.xi-riccatiPostEvent.chi);
//     stoPolicy.dts0   = - (1.0/sgm) * (riccatiPostEvent.eta-riccatiPostEvent.iota);
//     riccatiPreEvent.s.noalias()   
//         += (1.0/sgm) * (riccatiPostEvent.Psi-riccatiPostEvent.Phi) * (riccatiPostEvent.eta-riccatiPostEvent.iota);
//     riccatiPreEvent.Phi.noalias() 
//         -= (1.0/sgm) * (riccatiPostEvent.Psi-riccatiPostEvent.Phi) * (riccatiPostEvent.xi-riccatiPostEvent.chi);
//     riccatiPreEvent.rho
//         = riccatiPostEvent.xi - (1.0/sgm) * (riccatiPostEvent.xi-riccatiPostEvent.chi) * (riccatiPostEvent.xi-riccatiPostEvent.chi);
//     riccatiPreEvent.iota
//         = riccatiPostEvent.eta - (1.0/sgm) * (riccatiPostEvent.xi-riccatiPostEvent.chi)  * (riccatiPostEvent.eta-riccatiPostEvent.iota);
//   }
// }

} // namespace stoc
} // namespace ocs2