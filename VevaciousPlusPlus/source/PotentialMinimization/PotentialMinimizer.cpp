/*
 * PotentialMinimizer.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  PotentialMinimizer::PotentialMinimizer(
                                       PotentialFunction& potentialFunction ) :
    potentialFunction( potentialFunction ),
    foundMinima(),
    dsbVacuum(),
    panicVacua(),
    panicVacuum()
  {
    // This constructor is just an initialization list.
  }

  PotentialMinimizer::~PotentialMinimizer()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
