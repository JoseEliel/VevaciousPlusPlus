/*
 * HomotopyContinuationSolver.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/HomotopyContinuation/HomotopyContinuationSolver.hpp"

namespace VevaciousPlusPlus
{

  HomotopyContinuationSolver::HomotopyContinuationSolver(
      HomotopyContinuationTargetSystem const& homotopyContinuationPotential ) :
    StartingPointFinder(),
    homotopyContinuationPotential( homotopyContinuationPotential )
  {
    // This constructor is just an initialization list.
  }

  HomotopyContinuationSolver::~HomotopyContinuationSolver()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
