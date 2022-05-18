#pragma once

#include <memory>

#include <ocs2_sto_ipm/Types.h>
#include <ocs2_sto_ipm/interior_point_method/InteriorPointMethodData.h>

namespace ocs2 {
namespace sto_ipm {

void initInteriorPointMethodData(const VectorFunctionLinearApproximation& constraint,
                                 InteriorPointMethodData& data, scalar_t barrier);

void evalPerturbedResidual(const VectorFunctionLinearApproximation& constraint,
                           InteriorPointMethodData& data);

void evalLagrangianDerivativesStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateConstraint, 
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost);

void evalLagrangianDerivativesInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputConstraint, 
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost);

void evalLagrangianDerivativesStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputConstraint, 
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost);

void computeCondensingCoefficient(InteriorPointMethodData& data);

void condensingStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateIneqConstraint,
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost);

void condensingInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputIneqConstraint,
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost);

void condensingStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputIneqConstraint,
    InteriorPointMethodData& data, 
    ScalarFunctionQuadraticApproximation& cost);

void expansionStateIneqConstraint(
    const VectorFunctionLinearApproximation& stateIneqConstraint,
    const vector_t& dx, InteriorPointMethodData& data);

void expansionInputIneqConstraint(
    const VectorFunctionLinearApproximation& inputIneqConstraint,
    const vector_t& dx, InteriorPointMethodData& data);

void expansionStateInputIneqConstraint(
    const VectorFunctionLinearApproximation& stateInputIneqConstraint, 
    const vector_t& dx, const vector_t& du, InteriorPointMethodData& data);

scalar_t fractionToBoundaryStepSize(const InteriorPointMethodData& data);

}  // namespace sto_ipm
}  // namespace ocs2
