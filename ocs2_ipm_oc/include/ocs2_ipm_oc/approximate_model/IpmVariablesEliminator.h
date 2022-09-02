#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/ModelData.h>

#include "ocs2_ipm_oc/oc_data/IpmVariables.h"
#include "ocs2_ipm_oc/oc_data/IpmData.h"

namespace ocs2 {
namespace ipm {

/**
 * Eliminates the IPM-related variables and inequality constraints from a discrete-time LQ approximation.
 *
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in, out] modelData: The discrete-time LQ approximation model data.
 * @param [out] ipmData: The IPM-related data.
 * @param [out] barrierParam: The barrier parameter of the IPM.
 */
void eliminateIpmVariablesIntermediateLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData, scalar_t barrierParam);

/**
 * Eliminates the IPM-related variables and inequality constraints from a discrete-time LQ approximation.
 *
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in, out] modelData: The discrete-time LQ approximation model data.
 * @param [out] barrierParam: The barrier parameter of the IPM.
 * @return The output IPM data.
 */
IpmData eliminateIpmVariablesIntermediateLQ(const IpmVariables& ipmVariables, ModelData& modelData, scalar_t barrierParam);

/**
 * Eliminates the IPM-related variables and inequality constraints from a discrete-time LQ approximation at a jump event time.
 *
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in, out] modelData: The discrete-time LQ approximation model data.
 * @param [out] ipmData: The IPM-related data.
 * @param [out] barrierParam: The barrier parameter of the IPM.
 */
void eliminateIpmVariablesPreJumpLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData, scalar_t barrierParam);

/**
 * Eliminates the IPM-related variables and inequality constraints from a discrete-time LQ approximation at a jump event time.
 *
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in, out] modelData: The discrete-time LQ approximation model data.
 * @param [out] barrierParam: The barrier parameter of the IPM.
 * @return The output IPM data.
 */
IpmData eliminateIpmVariablesPreJumpLQ(const IpmVariables& ipmVariables, ModelData& modelData, scalar_t barrierParam);

/**
 * Eliminates the IPM-related variables and inequality constraints from a discrete-time LQ approximation at final time.
 *
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in, out] modelData: The discrete-time LQ approximation model data.
 * @param [out] ipmData: The IPM-related data.
 * @param [out] barrierParam: The barrier parameter of the IPM.
 */
void eliminateIpmVariablesFinalLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData, scalar_t barrierParam);

/**
 * Eliminates the IPM-related variables and inequality constraints from a discrete-time LQ approximation at final time.
 *
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in, out] modelData: The discrete-time LQ approximation model data.
 * @param [out] barrierParam: The barrier parameter of the IPM.
 * @return The output IPM data.
 */
IpmData eliminateIpmVariablesFinalLQ(const IpmVariables& ipmVariables, ModelData& modelData, scalar_t barrierParam);

}  // namespace ipm
}  // namespace ocs2
