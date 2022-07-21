#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/ModelData.h>

namespace ocs2 {
namespace ipm {

void projectIntermediateLQ(const ModelData& modelData, VectorFunctionLinearApproximation& constraintProjection, 
                           ModelData& projectedModelData);

}  // namespace ipm
}  // namespace ocs2