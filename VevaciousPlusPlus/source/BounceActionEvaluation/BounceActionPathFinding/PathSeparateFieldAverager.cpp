/*
 * PathSeparateFieldAverager.cpp
 *
 *  Created on: Sep 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/PathSeparateFieldAverager.hpp"

namespace VevaciousPlusPlus
{

  PathSeparateFieldAverager::PathSeparateFieldAverager(
                    BounceActionCalculator const* const bounceActionCalculator,
                                             unsigned int const minuitStrategy,
                                          double const minuitToleranceFraction,
                                           size_t const movesPerImprovement ) :
    MinuitBetweenPaths( bounceActionCalculator,
                            minuitStrategy,
                            minuitToleranceFraction,
                            movesPerImprovement )
  {
    // This constructor is just an initialization list.
  }

  PathSeparateFieldAverager::~PathSeparateFieldAverager()
  {
    // This does nothing.
  }


  // This takes minuitParameters as the set of fractions for each field which
  // should be taken as weightings for the straightPath node's value for that
  // field while the weighting for the field from the node from curvedPath is
  // 1.0 - that weighting, then it returns the path of straight lines between
  // the composed nodes.
  TunnelPath const* PathSeparateFieldAverager::PathForParameterization(
                          std::vector< double > const& minuitParameters ) const
  {
    std::vector< std::vector< double > > composedPath( straightPath );
    for( size_t nodeIndex( 1 );
         nodeIndex < numberOfSegments;
         ++nodeIndex )
    {
      std::vector< double > const& straightNode( straightPath[ nodeIndex ] );
      std::vector< double > const& curvedNode( (*curvedPath)[ nodeIndex ] );
      std::vector< double >& currentNode( composedPath[ nodeIndex ] );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        currentNode[ fieldIndex ]
        = ( ( minuitParameters[ fieldIndex ] * straightNode[ fieldIndex ] )
            + ( ( 1.0 - minuitParameters[ fieldIndex ] )
                * curvedNode[ fieldIndex ] ) );
      }
    }
    return new LinearSplineThroughNodes( composedPath,
                                         std::vector< double >( 0 ),
                                         pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
