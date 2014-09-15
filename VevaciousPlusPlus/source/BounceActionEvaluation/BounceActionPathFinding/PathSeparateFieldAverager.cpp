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
                          BounceActionCalculator* const bounceActionCalculator,
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
    std::vector< std::vector< double > > composedPath( *curvedPath );
    size_t const numberOfVaryingNodes( composedPath.size() - 2 );
    double const segmentFractionStep( 1.0
                         / static_cast< double >( numberOfVaryingNodes + 1 ) );
    std::vector< double > const& falseVacuum( composedPath.front() );
    std::vector< double > const& trueVacuum( composedPath.back() );
    for( size_t nodeIndex( 1 );
         nodeIndex <= numberOfVaryingNodes;
         ++nodeIndex )
    {
      double const trueStraightFraction( nodeIndex * segmentFractionStep );
      double const falseStraightFraction( 1.0 - trueStraightFraction );
      std::vector< double > const& curvedNode( (*curvedPath)[ nodeIndex ] );
      std::vector< double >& currentNode( composedPath[ nodeIndex ] );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        double const straightFieldValue(
                          ( falseStraightFraction * falseVacuum[ fieldIndex ] )
                       + ( trueStraightFraction * trueVacuum[ fieldIndex ] ) );
        currentNode[ fieldIndex ]
        = ( ( minuitParameters[ fieldIndex ] * straightFieldValue )
            + ( ( 1.0 - minuitParameters[ fieldIndex ] )
                * curvedNode[ fieldIndex ] ) );
      }
    }
    return new LinearSplineThroughNodes( composedPath,
                                         std::vector< double >( 0 ),
                                         pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
