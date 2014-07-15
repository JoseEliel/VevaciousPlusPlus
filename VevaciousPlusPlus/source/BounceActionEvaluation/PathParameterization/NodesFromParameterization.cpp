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
                                      std::vector< double > const& falseVacuum,
                                       std::vector< double > const& trueVacuum,
                                     size_t const numberOfIntermediateNodes ) :
    numberOfFields( falseVacuum.size() ),
    numberOfIntermediateNodes( numberOfIntermediateNodes ),
    pathNodes( numberOfIntermediateNodes,
               std::vector< double >( numberOfFields ) ),
    zeroParameterization( ( numberOfFields - 1 ),
                          0.0 ),
    initialStepSizes( zeroParameterization )
  {
    pathNodes.front() = falseVacuum;
    pathNodes.back() = trueVacuum;
    SetInitialParameterizationAndStepSizes( zeroParameterization,
                                            initialStepSizes );
  }

  NodesFromParameterization::~NodesFromParameterization()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
