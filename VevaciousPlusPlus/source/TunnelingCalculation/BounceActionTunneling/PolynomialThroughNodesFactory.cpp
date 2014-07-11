/*
 * PolynomialThroughNodesFactory.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/PolynomialThroughNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  PolynomialThroughNodesFactory::PolynomialThroughNodesFactory(
                                            std::string const& xmlArguments ) :
    PathFromNodesFactory( xmlArguments ),
    pathStepsInverse()
  {
    size_t const
    totalNumberOfNodes( nodesFromParameterization->NumberOfPathNodes() );
    double const stepSize( 1.0 / (double)( totalNumberOfNodes - 1 ) );
    Eigen::MatrixXd pathSteps( totalNumberOfNodes,
                               numberOfFields );
    for( size_t rowIndex( 0 );
         rowIndex < totalNumberOfNodes;
         ++rowIndex )
    {
      pathSteps( rowIndex,
                 0 ) = 1.0;
    }
    for( size_t columnIndex( 1 );
         columnIndex < numberOfFields;
         ++columnIndex )
    {
      pathSteps( 0,
                 columnIndex ) = 0.0;
    }
    for( size_t columnIndex( 1 );
         columnIndex < numberOfFields;
         ++columnIndex )
    {
      pathSteps( 0,
                 ( numberOfFields - 1 ) ) = 1.0;
    }
    for( size_t rowIndex( 1 );
         rowIndex < ( totalNumberOfNodes - 1 );
         ++rowIndex )
    {
      for( size_t columnIndex( 1 );
           columnIndex < numberOfFields;
           ++columnIndex )
      {
        pathSteps( rowIndex,
                   columnIndex ) = pow( ( (double)rowIndex * stepSize ),
                                        columnIndex );
      }
    }
    pathStepsInverse = pathSteps.inverse();
  }

  PolynomialThroughNodesFactory::~PolynomialThroughNodesFactory()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
