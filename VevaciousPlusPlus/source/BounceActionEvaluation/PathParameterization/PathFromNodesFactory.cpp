/*
 * PathFromNodesFactory.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/PathFromNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  PathFromNodesFactory::PathFromNodesFactory( size_t const numberOfFields,
                       NodesFromParameterization* nodesFromParameterization ) :
    TunnelPathFactory( numberOfFields,
                       std::vector< double >( numberOfFields,
                                              0.0 ) ),
    nodesFromParameterization( nodesFromParameterization )
  {
    // We need to update zeroParameterization to how nodesFromParameterization
    // parameterizes the "zero parameterization".
    nodesFromParameterization->SetInitialParameterizationAndStepSizes(
                                                          zeroParameterization,
                                                            initialStepSizes );
  }

  PathFromNodesFactory::~PathFromNodesFactory()
  {
    delete nodesFromParameterization;
  }

} /* namespace VevaciousPlusPlus */
