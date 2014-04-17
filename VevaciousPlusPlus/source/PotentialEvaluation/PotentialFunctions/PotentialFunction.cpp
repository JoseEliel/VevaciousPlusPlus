/*
 * PotentialFunction.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  PotentialFunction::PotentialFunction( SlhaManager& slhaManager ) :
    SlhaUpdatePropagator( slhaManager ),
    fieldNames(),
    numberOfFields( 0 ),
    dsbFieldValueInputs()
  {
    // This constructor is just an initialization list.
  }

  PotentialFunction::PotentialFunction( PotentialFunction const& copySource ) :
    SlhaUpdatePropagator( copySource.slhaManager ),
    fieldNames( copySource.fieldNames ),
    numberOfFields( copySource.numberOfFields ),
    dsbFieldValueInputs( copySource.dsbFieldValueInputs )
  {
    // This constructor is just an initialization list.
  }

  PotentialFunction::~PotentialFunction()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
