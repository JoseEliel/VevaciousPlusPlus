/*
 * FullPathVaryingMinuit.cpp
 *
 *  Created on: Jul 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/FullPathVaryingMinuit.hpp"

namespace VevaciousPlusPlus
{

  FullPathVaryingMinuit::FullPathVaryingMinuit( TunnelPathFactory* pathFactory,
                                              size_t const movesPerImprovement,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinuitPathFinder( movesPerImprovement,
                      minuitStrategy,
                      minuitToleranceFraction ),
    pathFactory( pathFactory ),
    currentMinuitResult()
  {
    // This constructor is just an initialization list.
  }

  FullPathVaryingMinuit::~FullPathVaryingMinuit()
  {
    delete pathFactory;
  }

} /* namespace VevaciousPlusPlus */
