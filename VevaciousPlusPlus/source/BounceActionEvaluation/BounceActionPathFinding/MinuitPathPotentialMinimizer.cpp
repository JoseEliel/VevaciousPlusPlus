/*
 * MinuitPathPotentialMinimizer.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinuitPathPotentialMinimizer.hpp"

namespace VevaciousPlusPlus
{

  MinuitPathPotentialMinimizer::MinuitPathPotentialMinimizer(
                                          TunnelPathFactory* const pathFactory,
                                    PotentialFunction const& potentialFunction,
                                    size_t const numberOfPotentialSamplePoints,
                                              size_t const movesPerImprovement,
                                             unsigned int const minuitStrategy,
                                      double const minuitToleranceFraction  ) :
    FullPathVaryingMinuit( pathFactory,
                           movesPerImprovement,
                           minuitStrategy,
                           minuitToleranceFraction ),
    potentialFunction( potentialFunction ),
    numberOfPotentialSamplePoints( numberOfPotentialSamplePoints ),
    numberOfFields( potentialFunction.NumberOfFieldVariables() ),
    pathSegmentSize( 1.0
                 / static_cast< double >( numberOfPotentialSamplePoints + 1 ) )
  {
    // This constructor is just an initialization list.
  }

  MinuitPathPotentialMinimizer::~MinuitPathPotentialMinimizer()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
