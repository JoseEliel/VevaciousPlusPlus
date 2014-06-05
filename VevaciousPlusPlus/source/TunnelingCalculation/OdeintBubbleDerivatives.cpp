/*
 * OdeintBubbleDerivatives.cpp
 *
 *  Created on: Jun 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{

  OdeintBubbleDerivatives::OdeintBubbleDerivatives(
                       PathFieldsAndPotential const& pathFieldsAndPotential ) :
    potentialDerivative(),
    firstDerivatives( pathFieldsAndPotential.FieldDerivatives() ),
    numberOfFields( pathFieldsAndPotential.FieldDerivatives().size() ),
    secondDerivatives( numberOfFields ),
    dampingFactor( pathFieldsAndPotential.DampingFactor() )
  {
    potentialDerivative.BecomeFirstDerivativeOf(
                             pathFieldsAndPotential.PotentialApproximation() );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < firstDerivatives.size();
         ++fieldIndex )
    {
      secondDerivatives[ fieldIndex ].BecomeFirstDerivativeOf(
                                              firstDerivatives[ fieldIndex ] );
    }
  }

  OdeintBubbleDerivatives::~OdeintBubbleDerivatives()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
