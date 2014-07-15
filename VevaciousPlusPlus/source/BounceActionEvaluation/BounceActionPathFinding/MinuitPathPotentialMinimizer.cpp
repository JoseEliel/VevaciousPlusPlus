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
                                                TunnelPathFactory* pathFactory,
                                    PotentialFunction const& potentialFunction,
                                           std::string const& xmlArguments  ) :
    FullPathVaryingMinuit( pathFactory,
                           xmlArguments ),
    potentialFunction( potentialFunction ),
    numberOfPotentialSamplePoints( 15 )
  {
    // This constructor is just an initialization list.
  }

  MinuitPathPotentialMinimizer::~MinuitPathPotentialMinimizer()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
