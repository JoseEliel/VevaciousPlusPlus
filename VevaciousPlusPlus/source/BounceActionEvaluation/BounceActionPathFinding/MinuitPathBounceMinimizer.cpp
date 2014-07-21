/*
 * MinuitPathBounceMinimizer.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinuitPathBounceMinimizer.hpp"

namespace VevaciousPlusPlus
{

  MinuitPathBounceMinimizer::MinuitPathBounceMinimizer(
                                          TunnelPathFactory* const pathFactory,
                          BounceActionCalculator* const bounceActionCalculator,
                                              size_t const movesPerImprovement,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    FullPathVaryingMinuit( pathFactory,
                           movesPerImprovement,
                           minuitStrategy,
                           minuitToleranceFraction ),
    bounceActionCalculator( bounceActionCalculator )
  {
    // This constructor is just an initialization list.
  }

  MinuitPathBounceMinimizer::~MinuitPathBounceMinimizer()
  {
    delete bounceActionCalculator;
  }

} /* namespace VevaciousPlusPlus */
