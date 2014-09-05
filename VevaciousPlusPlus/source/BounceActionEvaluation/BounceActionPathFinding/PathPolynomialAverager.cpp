/*
 * PathPolynomialAverager.cpp
 *
 *  Created on: Sep 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/PathPolynomialAverager.hpp"

namespace VevaciousPlusPlus
{

  PathPolynomialAverager::PathPolynomialAverager(
                                               size_t const degreeOfPolynomial,
                          BounceActionCalculator* const bounceActionCalculator,
                                             unsigned int const minuitStrategy,
                                          double const minuitToleranceFraction,
                                           size_t const movesPerImprovement ) :
    MinuitBetweenPaths( bounceActionCalculator,
                        minuitStrategy,
                        minuitToleranceFraction,
                        movesPerImprovement ),
    degreeOfPolynomial( degreeOfPolynomial )
  {
    // This constructor is just an initialization list.
  }

  PathPolynomialAverager::~PathPolynomialAverager()
  {
    // This does nothing.
  }


  // This takes minuitParameters as a set of coefficients for weightings and
  // creates a node path of weighted averages of each node in curvedPath with
  // the corresponding node in straightPath, where the weighting for the
  // straightPath node is
  // minuitParameters[ 0 ] + ( minuitParameters[ 1 ] * the fraction along the
  // path of the node) + ( minuitParameters[ 2 ] * that fraction squared )
  // + ..., and the weighting for the curvedPath node is 1.0 - the other
  // weighting, then it returns the path of straight lines between the
  // composed nodes.
  TunnelPath const* PathPolynomialAverager::PathForParameterization(
                          std::vector< double > const& minuitParameters ) const
  {
    std::vector< std::vector< double > > composedPath( straightPath );
    for( size_t nodeIndex( 1 );
         nodeIndex < numberOfSegments;
         ++nodeIndex )
    {
      std::vector< double > const& straightNode( straightPath[ nodeIndex ] );
      std::vector< double > const& curvedNode( (*curvedPath)[ nodeIndex ] );
      double const pathFraction( static_cast< double >( nodeIndex )
                                 / static_cast< double >( numberOfSegments ) );
      double straightFraction( minuitParameters[ 0 ] );
      for( size_t coefficientIndex( 1 );
           coefficientIndex <= degreeOfPolynomial;
           ++coefficientIndex )
      {
        straightFraction += ( minuitParameters[ coefficientIndex ]
                              * pow( pathFraction,
                                     coefficientIndex ) );
      }
      double const curvedFraction( 1.0 - straightFraction );
      std::vector< double >& currentNode( composedPath[ nodeIndex ] );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        currentNode[ fieldIndex ]
        = ( ( straightFraction * straightNode[ fieldIndex ] )
            + ( curvedFraction * curvedNode[ fieldIndex ] ) );
      }
    }
    return new LinearSplineThroughNodes( composedPath,
                                         std::vector< double >( 0 ),
                                         pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
