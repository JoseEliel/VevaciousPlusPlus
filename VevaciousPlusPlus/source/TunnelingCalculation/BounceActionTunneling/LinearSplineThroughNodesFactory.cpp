/*
 * LinearSplineThroughNodesFactory.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/LinearSplineThroughNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  LinearSplineThroughNodesFactory::LinearSplineThroughNodesFactory(
                                      std::vector< double > const& falseVacuum,
                                       std::vector< double > const& trueVacuum,
                                            std::string const& xmlArguments ) :
    PathFromNodesFactory( falseVacuum,
                          trueVacuum,
                          xmlArguments )
  {
    // This constructor is just an initialization list.
  }

  LinearSplineThroughNodesFactory::~LinearSplineThroughNodesFactory()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
