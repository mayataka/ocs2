#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/ModelData.h>

namespace ocs2 {
namespace ipm {

/**
 * Computes the constraint projection and projects the LQ model data along with the state-input equality constraint.
 * @param [in] modelData: The discrete-time LQ approximation model data.
 * @param [out] constraintProjection: Constraint projection
 * @param [out] projectedModelData: Output projected LQ model data.
 */
void projectIntermediateLQ(const ModelData& modelData, VectorFunctionLinearApproximation& constraintProjection, 
                           ModelData& projectedModelData);

}  // namespace ipm
}  // namespace ocs2