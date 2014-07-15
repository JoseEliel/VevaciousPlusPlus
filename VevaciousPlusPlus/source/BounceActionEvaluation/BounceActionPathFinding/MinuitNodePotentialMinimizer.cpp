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
                                    PotentialFunction const& potentialFunction,
                                             PathFromNodesFactory* pathFactory,
                                           std::string const& xmlArguments  ) :
    SingleNodeVaryingMinuit( potentialFunction,
                             pathFactory,
                             xmlArguments )
  {
    // This constructor is just an initialization list.
  }

  MinuitNodePotentialMinimizer::~MinuitNodePotentialMinimizer()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
