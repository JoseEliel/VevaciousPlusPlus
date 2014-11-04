/*
 * MinimizingPotentialOnParallelPlanes.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinimizingPotentialOnParallelPlanes.hpp"

namespace VevaciousPlusPlus
{
  MinimizingPotentialOnParallelPlanes::MinimizingPotentialOnParallelPlanes(
                                    PotentialFunction const& potentialFunction,
                                             size_t const numberOfPathSegments,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinimizingPotentialOnHypersurfaces( potentialFunction,
                                        numberOfPathSegments,
                                        minuitStrategy,
                                        minuitToleranceFraction ),
    notYetProvidedPath( true ),
    planeDifferenceFraction( 1.0
                             / static_cast< double > ( numberOfPathSegments ) )
  {
    // This constructor is just an initialization list.
  }

  MinimizingPotentialOnParallelPlanes::~MinimizingPotentialOnParallelPlanes()
  {
    // This does nothing.
  }



  // This minimizes the potential on a series of hyperplanes, all
  // perpendicular to the vector difference of the vacua. It tries to avoid
  // "zig-zag-iness" from Minuit2 coming to rest on a jittery line reasonably
  // far from the starting path by setting the starting point on each plane
  // to be the previous node plus the vector difference of previous node from
  // the node before it, with a special case for the first varying node. It
  // ignores both arguments, and also sets notYetProvidedPath to false.
  TunnelPath const* MinimizingPotentialOnParallelPlanes::TryToImprovePath(
                                                    TunnelPath const& lastPath,
                                      BubbleProfile const& bubbleFromLastPath )
  {
    SetParallelVector( pathNodes.front(),
                       pathNodes.back() );
    SetCurrentMinuitSteps( planeDifferenceFraction );
    SetUpHouseholderReflection();

    currentHyperplaneOrigin = pathNodes.front();
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentHyperplaneOrigin[ fieldIndex ] = ( pathNodes.front()[ fieldIndex ]
                                                + ( planeDifferenceFraction
                                  * currentParallelComponent[ fieldIndex ] ) );
    }
    RunMigradAndPutTransformedResultIn( pathNodes[ 1 ] );

    for( size_t nodeIndex( 2 );
         nodeIndex <= numberOfVaryingNodes;
         ++nodeIndex )
    {
      // Each minimization in the next hyperplane starts at the previous node
      // plus the difference vector going to the previous node from the node
      // previous to that
      // (e.g. node[ 3 ] = vector sum of node[ 2 ] + ( node[ 2 ] - node[ 1 ] ),
      // hence the factor of 2 and the minus sign).
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        currentHyperplaneOrigin[ fieldIndex ]
        = ( ( 2.0 * pathNodes[ nodeIndex - 1 ][ fieldIndex ] )
            - pathNodes[ nodeIndex - 2 ][ fieldIndex ] );
      }
      RunMigradAndPutTransformedResultIn( pathNodes[ nodeIndex ] );
    }
    notYetProvidedPath = false;
    return new LinearSplineThroughNodes( pathNodes,
                                         nodeZeroParameterization,
                                         pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
