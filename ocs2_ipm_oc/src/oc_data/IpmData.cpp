#include "ocs2_ipm_oc/oc_data/IpmData.h"

namespace ocs2 {
namespace ipm {

void initIpmData(const ModelData& modelData, IpmData& ipmData) {
  initInteriorPointMethodData(modelData.stateIneqConstraint, ipmData.dataStateIneqConstraint);
  initInteriorPointMethodData(modelData.stateInputIneqConstraint, ipmData.dataStateInputIneqConstraint);
}

IpmData initIpmData(const ModelData& modelData) {
  IpmData ipmData;
  initIpmData(modelData, ipmData);
  return ipmData;
}

}  // namespace ipm
}  // namespace ocs2
