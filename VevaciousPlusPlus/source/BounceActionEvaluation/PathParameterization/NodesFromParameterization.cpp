/*
 * NodesFromParameterization.cpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/NodesFromParameterization.hpp"

namespace VevaciousPlusPlus
{

  NodesFromParameterization::NodesFromParameterization(
                                                   size_t const numberOfFields,
                                     size_t const numberOfIntermediateNodes ) :
    numberOfFields( numberOfFields ),
    numberOfIntermediateNodes( numberOfIntermediateNodes ),
    pathNodes( ( numberOfIntermediateNodes + 2 ),
               std::vector< double >( numberOfFields ) ),
    zeroFullParameterization( ( ( numberOfFields - 1 )
                                * numberOfIntermediateNodes ),
                              0.0 ),
    zeroNodeParameterization( ( numberOfFields - 1 ),
                              0.0 ),
    initialStepSizes( zeroFullParameterization )
  {
    // This constructor is just an initialization list.
  }

  NodesFromParameterization::~NodesFromParameterization()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
