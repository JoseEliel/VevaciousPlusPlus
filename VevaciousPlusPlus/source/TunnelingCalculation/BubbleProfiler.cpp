/*
 * BubbleProfiler.cpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{

  BubbleProfiler::BubbleProfiler(
                                SimplePolynomial const& potentialApproximation,
                   std::vector< SimplePolynomial > const& fieldPathDerivatives,
                                  unsigned int const dampingFactor ) :
    potentialDerivative( potentialApproximation.FirstDerivative() ),
    numberOfFields( fieldPathDerivatives.size() ),
    firstDerivatives( fieldPathDerivatives ),
    secondDerivatives( numberOfFields ),
    dampingFactor( (double)dampingFactor )
  {
    for( unsigned int fieldIndex( 0 );
         fieldIndex < fieldPathDerivatives.size();
         ++fieldIndex )
    {
      secondDerivatives[ fieldIndex ]
      = fieldPathDerivatives[ fieldIndex ].FirstDerivative();
    }
  }

  BubbleProfiler::~BubbleProfiler()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
