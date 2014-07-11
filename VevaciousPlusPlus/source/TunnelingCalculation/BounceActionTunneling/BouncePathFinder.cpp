/*
 * BouncePathFinder.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/BouncePathFinder.hpp"

namespace VevaciousPlusPlus
{

  BouncePathFinder::BouncePathFinder() :
    pathTemperature( 0.0 ),
    pathCanBeImproved( true ),
    currentPath( NULL )
  {
    // This constructor is just an initialization list.
  }

  BouncePathFinder::~BouncePathFinder()
  {
    delete currentPath;
  }

} /* namespace VevaciousPlusPlus */
