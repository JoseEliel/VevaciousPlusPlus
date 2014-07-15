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
                                                TunnelPathFactory* pathFactory,
                                BounceActionCalculator* bounceActionCalculator,
                                           std::string const& xmlArguments  ) :
    FullPathVaryingMinuit( pathFactory,
                           xmlArguments ),
    bounceActionCalculator( bounceActionCalculator )
  {
    // This constructor is just an initialization list.
  }

  MinuitPathBounceMinimizer::~MinuitPathBounceMinimizer()
  {
    delete bounceActionCalculator;
  }

} /* namespace VevaciousPlusPlus */
