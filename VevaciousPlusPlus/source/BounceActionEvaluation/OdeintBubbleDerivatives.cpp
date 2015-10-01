/*
 * OdeintBubbleDerivatives.cpp
 *
 *  Created on: Jun 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/OdeintBubbleDerivatives.hpp"

namespace VevaciousPlusPlus
{

  OdeintBubbleDerivatives::OdeintBubbleDerivatives(
                         OneDimensionalPotentialAlongPath const& pathPotential,
                                               TunnelPath const& tunnelPath ) :
    pathPotential( pathPotential ),
    tunnelPath( tunnelPath ),
    dampingFactor( 3.0 )
  {
    if( tunnelPath.NonZeroTemperature() )
    {
      dampingFactor = 2.0;
    }
  }

  OdeintBubbleDerivatives::~OdeintBubbleDerivatives()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
