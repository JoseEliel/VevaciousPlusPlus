/*
 * HomotopyContinuationAndGradient.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  HomotopyContinuationAndGradient::HomotopyContinuationAndGradient(
                                    PotentialFunction const& potentialFunction,
                     HomotopyContinuationSolver& homotopyContinuationSolver ) :
    PotentialMinimizer( potentialFunction ),
    homotopyContinuationSolver( homotopyContinuationSolver ),
    purelyRealSolutionSets()
  {
    // This constructor is just an initialization list.
  }

  HomotopyContinuationAndGradient::~HomotopyContinuationAndGradient()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
