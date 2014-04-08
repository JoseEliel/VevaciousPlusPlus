/*
 * PotentialFunction.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  PotentialFunction::PotentialFunction( SlhaManager& slhaManager ) :
    BOL::BasicObserver(),
    fieldNames(),
    numberOfFields( 0 )
  {
    slhaManager.registerObserver( this );
  }

  PotentialFunction::~PotentialFunction()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
