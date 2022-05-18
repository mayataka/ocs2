#include <ocs2_sto_ipm/riccati_recursion/BackwardRiccatiRecursion.h>

#include <cmath>
#include <exception>

namespace ocs2 {
namespace sto_ipm {

BackwardRiccatiRecursion::BackwardRiccatiRecursion(const size_t nx, 
                                                   const size_t nu, 
                                                   const scalar_t dts0_max) 
  : nx_(nx),
    nu_(nu),
    dts0_max_(dts0_max),
    sgm_eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    AtP_(matrix_t::Zero(nx, nx)),
    BtP_(matrix_t::Zero(nu, nx)),
    GK_(matrix_t::Zero(nu, nx)),
    Ginv_3_(matrix_t::Zero(3, 3)),
    Ginv_2_(matrix_t::Zero(2, 2)),
    Ginv_1_(0),
    Pf_(vector_t::Zero(nx)),
    llt_(nu) {
  if (nx <= 0) {
    throw std::out_of_range("nx must be positive!");
  }
  if (nu < 0) {
    throw std::out_of_range("nu must be non-negative!");
  }
  if (dts0_max < 0) {
    throw std::out_of_range("dts0_max must be non-negative!");
  }
}


BackwardRiccatiRecursion::BackwardRiccatiRecursion() 
  : nx_(0),
    nu_(0),
    dts0_max_(0),
    sgm_eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    AtP_(),
    BtP_(),
    GK_(),
    Pf_(),
    Ginv_3_(),
    Ginv_2_(),
    Ginv_1_(0),
    llt_() {
}


void BackwardRiccatiRecursion::setRegularization(const scalar_t dts0_max) {
  if (dts0_max < 0) {
    throw std::out_of_range("dts0_max must be non-negative!");
  }
  dts0_max_ = dts0_max;
}


void BackwardRiccatiRecursion::compute(
    const ScalarFunctionQuadraticApproximationWrapper& cost,
    RiccatiRecursionData& riccati) {
  riccati.P = cost.base.dfdxx;
  riccati.s = - cost.base.dfdx;
}


void BackwardRiccatiRecursion::compute(
    const RiccatiRecursionData& riccati_next,
    const VectorFunctionLinearApproximationWrapper& dynamics, 
    ScalarFunctionQuadraticApproximationWrapper& cost, 
    RiccatiRecursionData& riccati, LQRPolicy& lqr_policy) {
  AtP_.noalias() = dynamics.base.dfdx.transpose() * riccati_next.P;
  BtP_.noalias() = dynamics.base.dfdu.transpose() * riccati_next.P;
  // Factorize F
  cost.base.dfdxx.noalias() += AtP_ * dynamics.base.dfdx;
  // Factorize H^T
  cost.base.dfdux.noalias() += BtP_ * dynamics.base.dfdx;
  // Factorize G
  cost.base.dfduu.noalias() += BtP_ * dynamics.base.dfdu;
  // Factorize vector term
  cost.base.dfdu.noalias() += BtP_ * dynamics.base.f;
  cost.base.dfdu.noalias() -= dynamics.base.dfdu.transpose() * riccati_next.s; 
  const int nu = cost.base.dfdu.size();
  switch (nu)
  {
  case 3:
    Ginv_3_ = cost.base.dfduu.inverse();
    lqr_policy.K.noalias() = - Ginv_3_ * cost.base.dfdux;
    lqr_policy.k.noalias() = - Ginv_3_ * cost.base.dfdu;
    break;
  case 2:
    Ginv_2_ = cost.base.dfduu.inverse();
    lqr_policy.K.noalias() = - Ginv_2_ * cost.base.dfdux;
    lqr_policy.k.noalias() = - Ginv_2_ * cost.base.dfdu;
    break;
  case 1:
    Ginv_1_ = 1.0 / cost.base.dfduu.coeff(0, 0);
    lqr_policy.K.noalias() = - Ginv_1_ * cost.base.dfdux;
    lqr_policy.k.noalias() = - Ginv_1_ * cost.base.dfdu;
    break;
  default:
    llt_.compute(cost.base.dfduu);
    lqr_policy.K.noalias() = - llt_.solve(cost.base.dfdux);
    lqr_policy.k.noalias() = - llt_.solve(cost.base.dfdu);
    break;
  }
}


void BackwardRiccatiRecursion::compute(
    const RiccatiRecursionData& riccati_next,
    const VectorFunctionLinearApproximationWrapper& dynamics, 
    ScalarFunctionQuadraticApproximationWrapper& cost, 
    RiccatiRecursionData& riccati, LQRPolicy& lqr_policy, 
    const bool sto, const bool has_next_sto_phase) {
  compute(riccati_next, dynamics, cost, riccati, lqr_policy);
  if (sto) {
    riccati.psi_x.noalias() = AtP_ * dynamics.dfdt;
    riccati.psi_u.noalias() = BtP_ * dynamics.dfdt;
    riccati.psi_x.noalias() += cost.dfdxt;
    riccati.psi_u.noalias() += cost.dfdut;
    riccati.psi_x.noalias() += dynamics.base.dfdx.transpose() * riccati_next.Psi;
    riccati.psi_u.noalias() += dynamics.base.dfdu.transpose() * riccati_next.Psi;
    if (has_next_sto_phase) {
      riccati.phi_x.noalias() = dynamics.base.dfdx.transpose() * riccati_next.Phi;
      riccati.phi_u.noalias() = dynamics.base.dfdu.transpose() * riccati_next.Phi;
    }
    else {
      riccati.phi_x.setZero();
      riccati.phi_u.setZero();
    }
    lqr_policy.T.noalias() = - llt_.solve(riccati.psi_u);
    if (has_next_sto_phase) {
      lqr_policy.W.noalias() = - llt_.solve(riccati.phi_u);
    }
    const int nu = cost.base.dfdu.size();
    switch (nu)
    {
    case 3:
      lqr_policy.T.noalias() = - Ginv_3_ * riccati.psi_u;
      if (has_next_sto_phase) {
        lqr_policy.W.noalias() = - Ginv_3_ * riccati.phi_u;
      }
      else {
        lqr_policy.W.setZero();
      }
      break;
    case 2:
      lqr_policy.T.noalias() = - Ginv_2_ * riccati.psi_u;
      if (has_next_sto_phase) {
        lqr_policy.W.noalias() = - Ginv_2_ * riccati.phi_u;
      }
      else {
        lqr_policy.W.setZero();
      }
      break;
    case 1:
      lqr_policy.T.noalias() = - Ginv_1_ * riccati.psi_u;
      if (has_next_sto_phase) {
        lqr_policy.W.noalias() = - Ginv_1_ * riccati.phi_u;
      }
      else {
        lqr_policy.W.setZero();
      }
      break;
    default:
      lqr_policy.T.noalias() = - llt_.solve(riccati.psi_u);
      if (has_next_sto_phase) {
        lqr_policy.W.noalias() = - llt_.solve(riccati.phi_u);
      }
      else {
        lqr_policy.W.setZero();
      }
      break;
    }
    // The cost-to-go w.r.t. switching times
    // Vtx
    riccati.Psi            = riccati.psi_x;
    riccati.Psi.noalias() += lqr_policy.K.transpose() * riccati.psi_u;
    if (has_next_sto_phase) {
      riccati.Phi.noalias()  = riccati.phi_x;
      riccati.Phi.noalias() += lqr_policy.K.transpose() * riccati.phi_u;
    }
    else {
      riccati.Phi.setZero();
    }
    // Vtt
    Pf_.noalias() = riccati_next.P * dynamics.dfdt;
    riccati.xi    = dynamics.dfdt.dot(Pf_);
    riccati.xi   += cost.dfdtt; 
    riccati.xi   += 2 * riccati_next.Psi.dot(dynamics.dfdt);
    riccati.xi   += lqr_policy.T.dot(riccati.psi_u);
    riccati.xi   += riccati_next.xi;
    if (has_next_sto_phase) {
      riccati.chi  = cost.dfdtt;
      // TODO: change this to 
      // riccati.chi  = cost.dfdtt_prev;
      riccati.chi += riccati_next.Phi.dot(dynamics.dfdt);
      riccati.chi += lqr_policy.T.dot(riccati.phi_u);
      riccati.chi += riccati_next.chi;
      riccati.rho  = lqr_policy.W.dot(riccati.phi_u);
      riccati.rho += riccati_next.rho;
    }
    else {
      riccati.chi = 0.0;
      riccati.rho = 0.0;
    }
    // Vt
    Pf_.noalias() = riccati_next.P * dynamics.base.f - riccati_next.s;
    riccati.eta   = dynamics.dfdt.dot(Pf_);
    riccati.eta  += cost.dfdt;
    riccati.eta  += riccati_next.Psi.dot(dynamics.base.f);
    riccati.eta  += riccati.psi_u.dot(lqr_policy.k);
    riccati.eta  += riccati_next.eta;
    if (has_next_sto_phase) {
      riccati.iota  = riccati_next.Phi.dot(dynamics.base.f);
      riccati.iota += riccati.phi_u.dot(lqr_policy.k);
      riccati.iota += riccati_next.iota;
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


void BackwardRiccatiRecursion::phaseTransition(
    const RiccatiRecursionData& riccati, 
    RiccatiRecursionData& riccati_pre_event, STOPolicy& sto_policy, 
    const bool has_next_sto_phase) const {
  riccati_pre_event.P = riccati.P;
  riccati_pre_event.s = riccati.s;
  riccati_pre_event.Psi.setZero();
  riccati_pre_event.Phi = riccati.Psi;
  riccati_pre_event.xi = 0.0;
  riccati_pre_event.chi = 0.0;
  riccati_pre_event.rho = riccati.xi;
  riccati_pre_event.eta = 0.0;
  riccati_pre_event.iota = riccati.eta;
  if (has_next_sto_phase) {
    double sgm = riccati.xi - 2.0 * riccati.chi + riccati.rho;
    if ((sgm*dts0_max_) < std::abs(riccati.eta-riccati.iota) || sgm < sgm_eps_) {
      sgm = std::abs(sgm) + std::abs(riccati.eta-riccati.iota) / dts0_max_;
    }
    sto_policy.dtsdx  = - (1.0/sgm) * (riccati.Psi-riccati.Phi);
    sto_policy.dtsdts =   (1.0/sgm) * (riccati.xi-riccati.chi);
    sto_policy.dts0   = - (1.0/sgm) * (riccati.eta-riccati.iota);
    riccati_pre_event.s.noalias()   
        += (1.0/sgm) * (riccati.Psi-riccati.Phi) * (riccati.eta-riccati.iota);
    riccati_pre_event.Phi.noalias() 
        -= (1.0/sgm) * (riccati.Psi-riccati.Phi) * (riccati.xi-riccati.chi);
    riccati_pre_event.rho
        = riccati.xi - (1.0/sgm) * (riccati.xi-riccati.chi) * (riccati.xi-riccati.chi);
    riccati_pre_event.iota
        = riccati.eta - (1.0/sgm) * (riccati.xi-riccati.chi)  * (riccati.eta-riccati.iota);
  }
}

} // namespace sto_ipm
} // namespace ocs2