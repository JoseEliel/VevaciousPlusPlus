/*
 * PolynomialPathFromCoefficients.cpp
 *
 *  Created on: Jul 21, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/PolynomialPathFromCoefficients.hpp"

namespace VevaciousPlusPlus
{

  PolynomialPathFromCoefficients::PolynomialPathFromCoefficients(
                                                   size_t const numberOfFields,
                                      std::vector< double > const& falseVacuum,
                                 std::vector< double > const& vacuumDifference,
                             std::vector< double > const& pathParameterization,
                                                  double const pathTemperature,
                                              size_t const referenceFieldIndex,
                                      size_t const coefficientsBeyondLinear ) :
    TunnelPath( numberOfFields,
                pathParameterization,
                pathTemperature )
  {
    std::vector< double > coefficientVector( 2 );
    coefficientVector[ 0 ] = falseVacuum[ referenceFieldIndex ];
    coefficientVector[ 1 ] = vacuumDifference[ referenceFieldIndex ];
    fieldPolynomials[ referenceFieldIndex ]
    = SimplePolynomial( coefficientVector );
    size_t const numberOfCoefficients( coefficientsBeyondLinear + 2 );
    coefficientVector.resize( numberOfCoefficients );
    for( size_t parameterizationIndex( 0 );
         parameterizationIndex < ( numberOfFields - 1 );
         ++parameterizationIndex )
    {
      size_t fieldIndex( parameterizationIndex );
      if( fieldIndex >= referenceFieldIndex )
      {
        ++fieldIndex;
      }
      coefficientVector[ 0 ] = falseVacuum[ fieldIndex ];
      coefficientVector[ 1 ] = vacuumDifference[ fieldIndex ];
      for( size_t coefficientIndex( 2 );
           coefficientIndex < numberOfCoefficients;
           ++coefficientIndex )
      {
        coefficientVector[ coefficientIndex ]
        = pathParameterization[ ( ( coefficientIndex - 2 ) * numberOfFields )
                                + parameterizationIndex ];
        coefficientVector[ 1 ] -= coefficientVector[ coefficientIndex ];
        // This adjusts the linear term so that the path comes back to the true
        // vacuum at p = 1.
      }
      fieldPolynomials[ fieldIndex ] = SimplePolynomial( coefficientVector );
    }
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      firstDerivatives[ fieldIndex ].BecomeFirstDerivativeOf(
                                              fieldPolynomials[ fieldIndex ] );
      secondDerivatives[ fieldIndex ].BecomeFirstDerivativeOf(
                                              firstDerivatives[ fieldIndex ] );
    }
  }

  PolynomialPathFromCoefficients::~PolynomialPathFromCoefficients()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
