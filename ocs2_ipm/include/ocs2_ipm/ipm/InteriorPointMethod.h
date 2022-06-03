#pragma once

#include <memory>

#include <ocs2_core/Types.h>
#include <ocs2_ipm/ipm/SlackDual.h>
#include <ocs2_ipm/ipm/InteriorPointMethodData.h>

namespace ocs2 {

void initSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                   SlackDual& slackDual, scalar_t barrier);

void initInteriorPointMethodData(const VectorFunctionLinearApproximation& ineqConstraint,
                                 InteriorPointMethodData& data);

void evalPerturbedResidual(const VectorFunctionLinearApproximation& ineqConstraint,
                           const SlackDual& slackDual, 
                           InteriorPointMethodData& data,
                           ScalarFunctionQuadraticApproximation& cost,
                           bool stateConstraint=true, bool inputConstraint=true);

void condenseSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                       const SlackDual& slackDual, InteriorPointMethodData& data, 
                       ScalarFunctionQuadraticApproximation& cost,
                       bool stateConstraint=true, bool inputConstraint=true);

void computeDualDirection(const SlackDual& slackDual, InteriorPointMethodData& data);

void expandSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                     const SlackDual& slackDual, const vector_t& dx, const vector_t& du,
                     InteriorPointMethodData& data);

void expandSlackDual(const VectorFunctionLinearApproximation& ineqConstraint,
                     const SlackDual& slackDual, const vector_t& dx, 
                     InteriorPointMethodData& data);

scalar_t fractionToBoundaryStepSize(size_t dim, const vector_t& v, const vector_t& dv,
                                    scalar_t marginRate=0.995);

scalar_t fractionToBoundaryPrimalStepSize(const SlackDual& slackDual, 
                                          const InteriorPointMethodData& data,
                                          scalar_t marginRate=0.995);

scalar_t fractionToBoundaryDualStepSize(const SlackDual& slackDual, 
                                        const InteriorPointMethodData& data,
                                        scalar_t marginRate=0.995);

scalar_t slackLogBarrier(const SlackDual& slackDual);

}  // namespace ocs2
