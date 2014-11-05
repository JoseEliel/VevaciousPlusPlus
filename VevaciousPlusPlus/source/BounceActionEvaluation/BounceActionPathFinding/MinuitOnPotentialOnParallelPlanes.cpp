/*
 * MinuitOnPotentialOnParallelPlanes.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialOnParallelPlanes.hpp"

namespace VevaciousPlusPlus
{
  MinuitOnPotentialOnParallelPlanes::MinuitOnPotentialOnParallelPlanes(
                                    PotentialFunction const& potentialFunction,
                                             size_t const numberOfPathSegments,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinuitOnHypersurfaces( potentialFunction,
                           numberOfPathSegments,
                           minuitStrategy,
                           minuitToleranceFraction )
  {
    // This constructor is just an initialization list.
  }

  MinuitOnPotentialOnParallelPlanes::~MinuitOnPotentialOnParallelPlanes()
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
  TunnelPath const* MinuitOnPotentialOnParallelPlanes::TryToImprovePath(
                                                    TunnelPath const& lastPath,
                                      BubbleProfile const& bubbleFromLastPath )
  {
    SetParallelVector( pathNodes.front(),
                       pathNodes.back() );
    SetCurrentMinuitSteps( segmentAuxiliaryLength );
    SetUpHouseholderReflection();

    currentHyperplaneOrigin = pathNodes.front();
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentHyperplaneOrigin[ fieldIndex ] = ( pathNodes.front()[ fieldIndex ]
                                                + ( segmentAuxiliaryLength
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
    return new LinearSplineThroughNodes( pathNodes,
                                         nodeZeroParameterization,
                                         pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
