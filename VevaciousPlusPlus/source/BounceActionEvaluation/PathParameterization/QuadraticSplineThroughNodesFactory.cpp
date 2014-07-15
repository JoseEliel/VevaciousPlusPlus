/*
 * QuadraticSplineThroughNodesFactory.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/QuadraticSplineThroughNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  QuadraticSplineThroughNodesFactory::QuadraticSplineThroughNodesFactory(
                                            std::string const& xmlArguments ) :
    PathFromNodesFactory( xmlArguments )
  {
    // This constructor is just an initialization list.
  }

  QuadraticSplineThroughNodesFactory::~QuadraticSplineThroughNodesFactory()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
