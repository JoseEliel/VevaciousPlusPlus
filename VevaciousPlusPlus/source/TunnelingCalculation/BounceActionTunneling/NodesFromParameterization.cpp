/*
 * NodesFromParameterization.cpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/NodesFromParameterization.hpp"

namespace VevaciousPlusPlus
{

  NodesFromParameterization::NodesFromParameterization(
                                      std::vector< double > const& falseVacuum,
                                       std::vector< double > const& trueVacuum,
                                     size_t const numberOfIntermediateNodes ) :
    numberOfFields( falseVacuum.size() ),
    numberOfIntermediateNodes( numberOfIntermediateNodes ),
    falseVacuum( falseVacuum ),
    trueVacuum( trueVacuum )
  {
    // This constructor is just an initialization list.
  }

  NodesFromParameterization::~NodesFromParameterization()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
