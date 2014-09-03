/*
 * MinuitBetweenPaths.cpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinuitBetweenPaths.hpp"

namespace VevaciousPlusPlus
{

  MinuitBetweenPaths::MinuitBetweenPaths(
                        std::vector< std::vector< double > > const& curvedPath,
                                          double const pathTemperature,
                 BounceActionCalculator const* const bounceActionCalculator ) :
    ROOT::Minuit2::FCNBase(),
    curvedPath( curvedPath ),
    pathTemperature( pathTemperature ),
    bounceActionCalculator( bounceActionCalculator ),
    straightPath(),
    numberOfFields( curvedPath.front().size() )
  {
    // This constructor is just an initialization list.
  }

  MinuitBetweenPaths::~MinuitBetweenPaths()
  {
    // This does nothing.
  }


  // This takes minuitParameters as a pair of weightings and creates a node
  // path of weighted averages of each node in curvedPath with the
  // corresponding node in straightPath, where the weighting for the
  // straightPath node is
  // minuitParameters[ 0 ] + ( minuitParameters[ 1 ] * the fraction along the
  // path of the node), and the weighting for the curvedPath node is
  // 1.0 - the other weighting, then it returns the path of straight lines
  // between the composed nodes.
  TunnelPath const* MinuitBetweenPaths::PathForParameterization(
                          std::vector< double > const& minuitParameters ) const
  {
    std::vector< std::vector< double > > composedPath( straightPath );
    for( size_t nodeIndex( 1 );
         nodeIndex < numberOfSegments;
         ++nodeIndex )
    {
      double const straightFraction( minuitParameters[ 0 ]
                                     + ( ( minuitParameters[ 1 ] * nodeIndex )
                               / static_cast< double >( numberOfSegments ) ) );
      double const curvedFraction( 1.0 - straightFraction );
      std::vector< double >& currentNode( composedPath[ nodeIndex ] );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        currentNode[ fieldIndex ]
        = ( ( straightFraction * straightPath[ fieldIndex ] )
            + ( curvedFraction * curvedPath[ fieldIndex ] ) );
      }
    }
    return new LinearSplineThroughNodes( composedPath,
                                         std::vector< double >( 0 ),
                                         pathTemperature );
  }

  // This sets up straightPath to go from curvedPath.front() to
  // curvedPath.back() in a straight line with as many nodes as curvedPath has,
  // as well as updating pathTemperature and currentMinuitTolerance.
  void MinuitBetweenPaths::UpdateNodes( double const pathTemperature )
  {
    this->pathTemperature = pathTemperature;
    straightPath.resize( curvedPath.size(),
                         curvedPath.front() );
    std::vector< double > straightNode( curvedPath.back() );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      straightNode[ fieldIndex ] -= curvedPath.front()[ fieldIndex ];
    }
    double const straightFraction( 1.0
                          / static_cast< double >( straightPath.size() - 1 ) );
    for( size_t nodeIndex( 1 );
         nodeIndex < ( straightPath.size() - 1 );
         ++nodeIndex )
    {
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        straightPath[ nodeIndex ][ fieldIndex ]
        += ( nodeIndex * straightFraction * straightNode[ fieldIndex ] );
      }
    }
    straightPath.back() = curvedPath.back();
    LinearSplineThroughNodes comparisonPath( curvedPath,
                                             std::vector< double >( 0 ),
                                             pathTemperature );
    currentMinuitTolerance = ( minuitToleranceFraction
                               * (*bounceActionCalculator)( comparisonPath ) );
  }

} /* namespace VevaciousPlusPlus */
