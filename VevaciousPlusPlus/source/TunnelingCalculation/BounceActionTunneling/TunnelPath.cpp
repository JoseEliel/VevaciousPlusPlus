/*
 * TunnelPath.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/TunnelPath.hpp"

namespace VevaciousPlusPlus
{

  TunnelPath::TunnelPath( size_t const numberOfFields,
                          double const temperatureValue ) :
    numberOfFields( numberOfFields ),
    temperatureValue( temperatureValue ),
    nonZeroTemperature( temperatureValue > 0.0 )
  {
    // This constructor is just an initialization list.
  }

  TunnelPath::~TunnelPath()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
