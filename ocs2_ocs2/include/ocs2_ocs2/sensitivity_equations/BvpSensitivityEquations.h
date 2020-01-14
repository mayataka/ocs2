/******************************************************************************
Copyright (c) 2017, Farbod Farshidian. All rights reserved.

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

#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <vector>

#include <ocs2_core/Dimensions.h>
#include <ocs2_core/integration/OdeBase.h>
#include <ocs2_core/misc/LinearInterpolation.h>
#include <ocs2_core/model_data/ModelDataLinearInterpolation.h>

namespace ocs2 {

/**
 * BVP sensitivity equations.
 *
 * @tparam STATE_DIM: Dimension of the state space.
 * @tparam INPUT_DIM: Dimension of the control input space.
 */
template <size_t STATE_DIM, size_t INPUT_DIM>
class BvpSensitivityEquations final : public OdeBase<STATE_DIM> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  using BASE = OdeBase<STATE_DIM>;

  using DIMENSIONS = Dimensions<STATE_DIM, INPUT_DIM>;
  using scalar_t = typename DIMENSIONS::scalar_t;
  using scalar_array_t = typename DIMENSIONS::scalar_array_t;
  using state_vector_t = typename DIMENSIONS::state_vector_t;
  using state_vector_array_t = typename DIMENSIONS::state_vector_array_t;
  using input_vector_t = typename DIMENSIONS::input_vector_t;
  using input_vector_array_t = typename DIMENSIONS::input_vector_array_t;
  using input_state_matrix_t = typename DIMENSIONS::input_state_matrix_t;
  using input_state_matrix_array_t = typename DIMENSIONS::input_state_matrix_array_t;
  using state_matrix_t = typename DIMENSIONS::state_matrix_t;
  using state_matrix_array_t = typename DIMENSIONS::state_matrix_array_t;
  using input_matrix_t = typename DIMENSIONS::input_matrix_t;
  using input_matrix_array_t = typename DIMENSIONS::input_matrix_array_t;
  using state_input_matrix_t = typename DIMENSIONS::state_input_matrix_t;
  using state_input_matrix_array_t = typename DIMENSIONS::state_input_matrix_array_t;
  using constraint1_vector_t = typename DIMENSIONS::constraint1_vector_t;
  using constraint1_vector_array_t = typename DIMENSIONS::constraint1_vector_array_t;
  using constraint1_state_matrix_t = typename DIMENSIONS::constraint1_state_matrix_t;
  using constraint1_state_matrix_array_t = typename DIMENSIONS::constraint1_state_matrix_array_t;
  using dynamic_vector_t = typename DIMENSIONS::dynamic_vector_t;
  using dynamic_matrix_t = typename DIMENSIONS::dynamic_matrix_t;

  /**
   * Constructor.
   */
  BvpSensitivityEquations() = default;

  /**
   * Default destructor.
   */
  ~BvpSensitivityEquations() override = default;

  /**
   * Returns pointer to the class.
   *
   * @return A raw pointer to the class.
   */
  BvpSensitivityEquations<STATE_DIM, INPUT_DIM>* clone() const { return new BvpSensitivityEquations<STATE_DIM, INPUT_DIM>(*this); }

  /**
   * Sets Data
   */
  void setData(const scalar_array_t* timeStampPtr, const ModelDataBase::array_t* modelDataPtr, const state_matrix_array_t* AmConstrainedPtr,
               const input_state_matrix_array_t* CmProjectedPtr, const state_vector_array_t* costatePtr,
               const constraint1_vector_array_t* lagrangianPtr, const scalar_array_t* controllerTimeStampPtr,
               const input_state_matrix_array_t* KmConstrainedPtr, const state_matrix_array_t* SmPtr) {
    timeStampPtr_ = timeStampPtr;
    modelDataPtr_ = modelDataPtr;
    AmConstrainedPtr_ = AmConstrainedPtr;
    CmProjectedPtr_ = CmProjectedPtr;
    costatePtr_ = costatePtr;
    lagrangianPtr_ = lagrangianPtr;
    controllerTimeStampPtr_ = controllerTimeStampPtr;
    KmConstrainedPtr_ = KmConstrainedPtr;
    SmPtr_ = SmPtr;
  }

  /**
   * Sets the multiplier of exogenous part of the equation. It is either zero
   * or plus-minus 1/(s_{i+1}-s_{i})
   *
   * @param [in] multiplier: the multiplier of exogenous part of the equation.
   */
  void setMultiplier(const scalar_t& multiplier) { multiplier_ = multiplier; }

  /**
   * Computes Derivative
   * @param [in] time: Normalized transition time
   * @param [in] Mv: transition state
   * @param [out] dMv: mapped state after transition
   */
  void computeFlowMap(const scalar_t& z, const state_vector_t& Mv, state_vector_t& dMvdz) override {
    BASE::numFunctionCalls_++;

    // denormalized time
    const scalar_t t = -z;
    auto indexAlpha = LinearInterpolation<scalar_t>::timeSegment(t, timeStampPtr_);

    ModelData::LinearInterpolation::interpolate(indexAlpha, Fv_, modelDataPtr_, ModelData::dynamics);
    ModelData::LinearInterpolation::interpolate(indexAlpha, Am_, modelDataPtr_, ModelData::dynamicsStateDerivative);
    ModelData::LinearInterpolation::interpolate(indexAlpha, Bm_, modelDataPtr_, ModelData::dynamicsInputDerivative);
    ModelData::LinearInterpolation::interpolate(indexAlpha, Qv_, modelDataPtr_, ModelData::costStateDerivative);
    ModelData::LinearInterpolation::interpolate(indexAlpha, Cm_, modelDataPtr_, ModelData::stateInputEqConstrStateDerivative);
    auto nc1 = (*modelDataPtr_)[indexAlpha.first].numStateInputEqConstr_;

    EigenLinearInterpolation<state_matrix_t>::interpolate(indexAlpha, AmConstrained_, AmConstrainedPtr_);
    EigenLinearInterpolation<input_state_matrix_t>::interpolate(indexAlpha, CmProjected_, CmProjectedPtr_);
    EigenLinearInterpolation<state_vector_t>::interpolate(indexAlpha, costate_, costatePtr_);
    EigenLinearInterpolation<constraint1_vector_t>::interpolate(indexAlpha, lagrangian_, lagrangianPtr_);

    indexAlpha = LinearInterpolation<scalar_t>::timeSegment(t, controllerTimeStampPtr_);
    EigenLinearInterpolation<input_state_matrix_t>::interpolate(indexAlpha, KmConstrained_, KmConstrainedPtr_);
    EigenLinearInterpolation<state_matrix_t>::interpolate(indexAlpha, Sm_, SmPtr_);

    // here we have used RmConstrained = (I-DmConstrained).transpose() * Rm
    // and Km = -(I-DmConstrained) \tilde{L} - CmProjected_
    dMvdz = (AmConstrained_ + Bm_ * (CmProjected_ + KmConstrained_)).transpose() * Mv +
            multiplier_ * (Qv_ + Am_.transpose() * costate_ + Cm_.topRows(nc1).transpose() * lagrangian_.head(nc1) + Sm_ * Fv_);
  }

 private:
  scalar_t multiplier_ = 0.0;

  const scalar_array_t* timeStampPtr_;
  const ModelDataBase::array_t* modelDataPtr_;
  const state_matrix_array_t* AmConstrainedPtr_;
  const input_state_matrix_array_t* CmProjectedPtr_;
  const state_vector_array_t* costatePtr_;
  const constraint1_vector_array_t* lagrangianPtr_;
  const scalar_array_t* controllerTimeStampPtr_;
  const input_state_matrix_array_t* KmConstrainedPtr_;
  const state_matrix_array_t* SmPtr_;

  dynamic_vector_t Fv_;
  dynamic_matrix_t Am_;
  dynamic_matrix_t Bm_;
  dynamic_vector_t Qv_;
  dynamic_matrix_t Cm_;
  state_matrix_t AmConstrained_;
  input_state_matrix_t CmProjected_;
  state_vector_t costate_;
  constraint1_vector_t lagrangian_;
  input_state_matrix_t KmConstrained_;
  state_matrix_t Sm_;
};

}  // namespace ocs2
