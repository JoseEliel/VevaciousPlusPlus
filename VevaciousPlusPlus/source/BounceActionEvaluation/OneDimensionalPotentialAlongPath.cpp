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
    auxiliaryOfPathFalseVacuum( 0.0 ),
    auxiliaryOfPathPanicVacuum( 1.0 ),
    definiteUndershootAuxiliary( 0.0 ),
    thresholdForNearPathPanic( 0.01 )
  {
    // This constructor is just an initialization list.
  }

  OneDimensionalPotentialAlongPath::~OneDimensionalPotentialAlongPath()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
