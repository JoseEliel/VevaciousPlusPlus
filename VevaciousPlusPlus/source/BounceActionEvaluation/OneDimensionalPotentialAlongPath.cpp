/*
 * OneDimensionalPotentialAlongPath.cpp
 *
 *  Created on: Nov 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/OneDimensionalPotentialAlongPath.hpp"

namespace VevaciousPlusPlus
{

  OneDimensionalPotentialAlongPath::OneDimensionalPotentialAlongPath() :
    definiteUndershootAuxiliary( -1.0 ),
    definiteOvershootAuxiliary( -1.0 ),
    thresholdForNearPathPanic( -1.0 )
  {
    // This constructor is just an initialization list.
  }

  OneDimensionalPotentialAlongPath::~OneDimensionalPotentialAlongPath()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
