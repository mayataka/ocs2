#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/core/InteriorPointMethodData.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>
#include <ocs2_ipm/model_data/ModelData.h>

namespace ocs2 {
namespace ipm {

/**
 * OC datas for the interior point method.
 */
struct IpmData {
  // Interior point method related datas
  InteriorPointMethodData dataStateIneqConstraint; 
  InteriorPointMethodData dataStateInputIneqConstraint; 
};

void initIpmData(const ModelData& modelData, IpmData& ipmData);

IpmData initIpmData(const ModelData& modelData);

}  // namespace ipm
}  // namespace ocs2
