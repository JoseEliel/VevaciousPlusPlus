/*
 * QuadraticSplineThroughNodesFactory.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/QuadraticSplineThroughNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  QuadraticSplineThroughNodesFactory::QuadraticSplineThroughNodesFactory(
                                      std::vector< double > const& falseVacuum,
                                       std::vector< double > const& trueVacuum,
                                            std::string const& xmlArguments ) :
    PathFromNodesFactory( falseVacuum,
                          trueVacuum,
                          xmlArguments )
  {
    // This constructor is just an initialization list.
  }

  QuadraticSplineThroughNodesFactory::~QuadraticSplineThroughNodesFactory()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
