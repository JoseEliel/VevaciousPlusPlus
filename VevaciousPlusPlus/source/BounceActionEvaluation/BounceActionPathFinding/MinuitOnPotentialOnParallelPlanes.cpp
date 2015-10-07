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
    SetParallelVector( returnPathNodes.front(),
                       returnPathNodes.back() );

    // A tenth of the Euclidean distance between hyperplanes is used as the
    // initial step size for Minuit, though the step size should be adapted
    // internally as the minimization proceeds.
    SetCurrentMinuitSteps( 0.1 * segmentAuxiliaryLength );
    SetUpHouseholderReflection();

    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentHyperplaneOrigin( fieldIndex )
      = ( returnPathNodes.front()[ fieldIndex ]
          + ( segmentAuxiliaryLength
              * currentParallelComponent( fieldIndex ) ) );
    }
    RunMigradAndPutTransformedResultIn( returnPathNodes[ 1 ] );

    for( size_t nodeIndex( 2 );
         nodeIndex <= numberOfVaryingNodes;
         ++nodeIndex )
    {
      // Each minimization in the next hyperplane starts at points along the
      // straight vector between the vacua at equal intervals.
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        currentHyperplaneOrigin( fieldIndex )
        += ( segmentAuxiliaryLength * currentParallelComponent( fieldIndex ) );
      }
      RunMigradAndPutTransformedResultIn( returnPathNodes[ nodeIndex ] );
    }
    return new LinearSplineThroughNodes( returnPathNodes,
                                         nodeZeroParameterization,
                                         pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
