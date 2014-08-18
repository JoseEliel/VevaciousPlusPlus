/*
 * PolynomialThroughNodesFactory.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/PolynomialThroughNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  PolynomialThroughNodesFactory::PolynomialThroughNodesFactory(
                                                   size_t const numberOfFields,
                 NodesFromParameterization* const nodesFromParameterization ) :
    PathFromNodesFactory( numberOfFields,
                          nodesFromParameterization ),
    pathStepsInverse()
  {
    size_t const
    totalNumberOfNodes( nodesFromParameterization->NumberOfPathNodes() );
    double const stepSize( 1.0
                           / static_cast< double >( totalNumberOfNodes - 1 ) );
    Eigen::MatrixXd pathSteps( totalNumberOfNodes,
                               totalNumberOfNodes );
    for( size_t rowIndex( 0 );
         rowIndex < totalNumberOfNodes;
         ++rowIndex )
    {
      pathSteps( rowIndex,
                 0 ) = 1.0;
    }
    for( size_t columnIndex( 1 );
         columnIndex < totalNumberOfNodes;
         ++columnIndex )
    {
      pathSteps( 0,
                 columnIndex ) = 0.0;
    }
    for( size_t columnIndex( 1 );
         columnIndex < totalNumberOfNodes;
         ++columnIndex )
    {
      pathSteps( ( totalNumberOfNodes - 1 ),
                 columnIndex ) = 1.0;
    }
    for( size_t rowIndex( 1 );
         rowIndex < ( totalNumberOfNodes - 1 );
         ++rowIndex )
    {
      for( size_t columnIndex( 1 );
           columnIndex < totalNumberOfNodes;
           ++columnIndex )
      {
        pathSteps( rowIndex,
                   columnIndex ) = pow( ( rowIndex * stepSize ),
                                        columnIndex );
      }
    }
    pathStepsInverse = pathSteps.inverse();

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "pathSteps = " << std::endl << pathSteps << std::endl
    << "pathStepsInverse = " << std::endl << pathStepsInverse;
    std::cout << std::endl;*/
  }

  PolynomialThroughNodesFactory::~PolynomialThroughNodesFactory()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
