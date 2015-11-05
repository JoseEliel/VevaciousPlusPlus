/*
 * RgeImprovedOneLoopPotential.cpp
 *
 *  Created on: Nov 5, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/PotentialFunctions/RgeImprovedOneLoopPotential.hpp"

namespace VevaciousPlusPlus
{

  RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential(
                                              std::string const& modelFilename,
                               double const assumedPositiveOrNegativeTolerance,
                       ParameterUpdatePropagator& parameterUpdatePropagator ) :
    PotentialFromPolynomialWithMasses( modelFilename,
                                       assumedPositiveOrNegativeTolerance,
                   parameterUpdatePropagator.GetLagrangianParameterManager() ),
    ParameterUpdatePropagator( parameterUpdatePropagator ),
    minimumScaleSquared( -1.0 )
  {
    // This constructor is just an initialization list.
  }

  RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential(
                      PotentialFromPolynomialWithMasses const& potentialToCopy,
                       ParameterUpdatePropagator& parameterUpdatePropagator ) :
    PotentialFromPolynomialWithMasses( potentialToCopy ),
    ParameterUpdatePropagator( parameterUpdatePropagator ),
    minimumScaleSquared( -1.0 )
  {
    // This constructor is just an initialization list.
  }

  RgeImprovedOneLoopPotential::~RgeImprovedOneLoopPotential()
  {
    // This does nothing.
  }


  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  double RgeImprovedOneLoopPotential::operator()(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
  {
    double scaleSquared( temperatureValue * temperatureValue );
    for( std::vector< double >::const_iterator
         fieldValue( fieldConfiguration.begin() );
         fieldValue < fieldConfiguration.end();
         ++fieldValue )
    {
      scaleSquared += ( (*fieldValue) * (*fieldValue) );
    }
    scaleSquared = std::min( scaleSquared, minimumScaleSquared );
    // The logarithm of the scale is of course half the logarithm of the square
    // of the scale.
    std::vector< double > const
    parameterValues( lagrangianParameterManager.ParameterValues(
                                             ( 0.5 * log( scaleSquared ) ) ) );
    std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( parameterValues,
                                      fieldConfiguration,
                                      scalarSquareMasses,
                                      scalarMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( parameterValues,
                                      fieldConfiguration,
                                      fermionSquareMasses,
                                      fermionMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( parameterValues,
                                      fieldConfiguration,
                                      vectorSquareMasses,
                                      vectorMassesSquaredWithFactors );
    return ( treeLevelPotential( parameterValues,
                                 fieldConfiguration )
             + polynomialLoopCorrections( parameterValues,
                                          fieldConfiguration )
             + LoopAndThermalCorrections( scalarMassesSquaredWithFactors,
                                          fermionMassesSquaredWithFactors,
                                          vectorMassesSquaredWithFactors,
                                          ( 1.0 / scaleSquared ),
                                          temperatureValue ) );
  }

} /* namespace VevaciousPlusPlus */
