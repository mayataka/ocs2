#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/core/InteriorPointMethodData.h>

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

}  // namespace ipm
}  // namespace ocs2
