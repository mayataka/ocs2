#pragma once

#include <memory>

#include <ocs2_core/Types.h>
#include <ocs2_ipm/ipm/SlackDualData.h>
#include <ocs2_ipm/ipm/InteriorPointMethodData.h>

namespace ocs2 {

void initSlackDualData(const VectorFunctionLinearApproximation& constraint,
                       SlackDualData& slackDual, scalar_t barrier);

void initInteriorPointMethodData(const VectorFunctionLinearApproximation& constraint,
                                 InteriorPointMethodData& data);

void evalPerturbedResidual(const VectorFunctionLinearApproximation& constraint,
                           const SlackDualData& slackDual, 
                           InteriorPointMethodData& data);

void evalLagrangianDerivativesStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateConstraint, 
    const SlackDualData& slackDual, ScalarFunctionQuadraticApproximation& cost);

void evalLagrangianDerivativesInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputConstraint, 
    const SlackDualData& slackDual, ScalarFunctionQuadraticApproximation& cost);

void evalLagrangianDerivativesStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputConstraint, 
    const SlackDualData& slackDual, ScalarFunctionQuadraticApproximation& cost);

void computeCondensingCoefficient(const SlackDualData& slackDual, InteriorPointMethodData& data);

void condensingStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateIneqConstraint,
    const SlackDualData& slackDual, InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost);

void condensingInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputIneqConstraint,
    const SlackDualData& slackDual, InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost);

void condensingStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputIneqConstraint,
    const SlackDualData& slackDual, InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost);

void computeDualDirection(const SlackDualData& slackDual, InteriorPointMethodData& data);

void expansionStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateIneqConstraint,
    const SlackDualData& slackDual, const vector_t& dx, InteriorPointMethodData& data);

void expansionInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputIneqConstraint,
    const SlackDualData& slackDual, const vector_t& du, InteriorPointMethodData& data);

void expansionStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputIneqConstraint, 
    const SlackDualData& slackDual, const vector_t& dx, const vector_t& du, 
    InteriorPointMethodData& data);

scalar_t fractionToBoundaryStepSize(size_t dim, const vector_t& v, const vector_t& dv,
                                    scalar_t marginRate=0.995);

scalar_t fractionToBoundaryPrimalStepSize(const SlackDualData& slackDual, 
                                          const InteriorPointMethodData& data,
                                          scalar_t marginRate=0.995);

scalar_t fractionToBoundaryDualStepSize(const SlackDualData& slackDual, 
                                        const InteriorPointMethodData& data,
                                        scalar_t marginRate=0.995);

}  // namespace ocs2
