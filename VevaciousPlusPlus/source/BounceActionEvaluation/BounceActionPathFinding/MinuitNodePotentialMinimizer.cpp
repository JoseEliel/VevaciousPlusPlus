/*
 * MinuitNodePotentialMinimizer.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinuitNodePotentialMinimizer.hpp"

namespace VevaciousPlusPlus
{

  MinuitNodePotentialMinimizer::MinuitNodePotentialMinimizer(
                                       PathFromNodesFactory* const pathFactory,
                                    PotentialFunction const& potentialFunction,
                                              size_t const movesPerImprovement,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinuitPathFinder( movesPerImprovement,
                      minuitStrategy,
                      minuitToleranceFraction )
  {
    // This constructor is just an initialization list.
  }

  MinuitNodePotentialMinimizer::~MinuitNodePotentialMinimizer()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
