#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/core/SlackDual.h>
#include <ocs2_ipm/core/SlackDualDirection.h>
#include <ocs2_ipm/core/InteriorPointMethodData.h>

namespace ocs2 {
namespace ipm {

void initSlackDual(const VectorFunctionLinearApproximation& ineqConstraint, SlackDual& slackDual, scalar_t barrierParam);

SlackDual initSlackDual(const VectorFunctionLinearApproximation& ineqConstraint, scalar_t barrierParam);

void initInteriorPointMethodData(const VectorFunctionLinearApproximation& ineqConstraint, InteriorPointMethodData& ipmData);

InteriorPointMethodData initInteriorPointMethodData(const VectorFunctionLinearApproximation& ineqConstraint);

void evalPerturbedResidual(const VectorFunctionLinearApproximation& ineqConstraint, const SlackDual& slackDual, 
                           InteriorPointMethodData& ipmData, ScalarFunctionQuadraticApproximation& cost, scalar_t barrierParam,
                           bool stateConstraint=true, bool inputConstraint=true);

void condenseSlackDual(const VectorFunctionLinearApproximation& ineqConstraint, const SlackDual& slackDual, 
                       InteriorPointMethodData& ipmData, ScalarFunctionQuadraticApproximation& cost,
                       bool stateConstraint=true, bool inputConstraint=true);

void expandDual(const SlackDual& slackDual, const InteriorPointMethodData& ipmData, 
                SlackDualDirection& slackDualDirection);

void expandSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                     const SlackDual& slackDual, const InteriorPointMethodData& ipmData, 
                     const vector_t& dx, const vector_t& du, SlackDualDirection& slackDualDirection);

void expandSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                     const SlackDual& slackDual, const InteriorPointMethodData& ipmData, 
                     const vector_t& dx, SlackDualDirection& slackDualDirection);

scalar_t fractionToBoundaryStepSize(size_t dim, const vector_t& v, const vector_t& dv,
                                    scalar_t marginRate=0.995);

scalar_t fractionToBoundaryPrimalStepSize(const SlackDual& slackDual, const SlackDualDirection& slackDualDirection,
                                          const InteriorPointMethodData& ipmData, scalar_t marginRate=0.995);

scalar_t fractionToBoundaryDualStepSize(const SlackDual& slackDual, const SlackDualDirection& slackDualDirection,
                                        const InteriorPointMethodData& ipmData, scalar_t marginRate=0.995);

}  // namespace ipm
}  // namespace ocs2
