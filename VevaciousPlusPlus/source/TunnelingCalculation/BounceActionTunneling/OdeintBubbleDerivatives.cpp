/*
 * OdeintBubbleDerivatives.cpp
 *
 *  Created on: Jun 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/OdeintBubbleDerivatives.hpp"

namespace VevaciousPlusPlus
{

  OdeintBubbleDerivatives::OdeintBubbleDerivatives(
                                          SplinePotential const& pathPotential,
                                               TunnelPath const& tunnelPath ) :
    pathPotential( pathPotential ),
    tunnelPath( tunnelPath ),
    dampingFactor( 3 )
  {
    if( tunnelPath.NonZeroTemperature() )
    {
      dampingFactor = 2;
    }
  }

  OdeintBubbleDerivatives::~OdeintBubbleDerivatives()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
