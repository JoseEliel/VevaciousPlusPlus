/*
 * NodesOnPlanes.cpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/NodesOnPlanes.hpp"

namespace VevaciousPlusPlus
{

  NodesOnPlanes::NodesOnPlanes( size_t const numberOfFields,
                                size_t const numberOfIntermediateNodes ) :
    NodesFromParameterization( numberOfFields,
                               numberOfIntermediateNodes ),
    numberOfParametersPerNode( numberOfFields - 1 )
  {
    // This constructor is just an initialization list.
  }

  NodesOnPlanes::~NodesOnPlanes()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
