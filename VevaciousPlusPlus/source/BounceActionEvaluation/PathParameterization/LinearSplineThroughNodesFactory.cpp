/*
 * LinearSplineThroughNodesFactory.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/LinearSplineThroughNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  LinearSplineThroughNodesFactory::LinearSplineThroughNodesFactory(
                                            std::string const& xmlArguments ) :
    PathFromNodesFactory( xmlArguments )
  {
    // This constructor is just an initialization list.
  }

  LinearSplineThroughNodesFactory::~LinearSplineThroughNodesFactory()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
