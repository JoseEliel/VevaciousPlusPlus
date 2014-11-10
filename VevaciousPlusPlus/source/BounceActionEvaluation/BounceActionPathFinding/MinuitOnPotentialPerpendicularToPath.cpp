/*
 * MinuitOnPotentialPerpendicularToPath.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialPerpendicularToPath.hpp"

namespace VevaciousPlusPlus
{

  MinuitOnPotentialPerpendicularToPath::MinuitOnPotentialPerpendicularToPath(
                                    PotentialFunction const& potentialFunction,
                                             size_t const numberOfPathSegments,
                                           int const numberOfAllowedWorsenings,
                                    double const nodeMovementThresholdFraction,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
     MinuitOnHypersurfaces( potentialFunction,
                            numberOfPathSegments,
                            minuitStrategy,
                            minuitToleranceFraction ),
     numberOfAllowedWorsenings( numberOfAllowedWorsenings ),
     nodeMovementThresholdFractionSquared( 0.25 * nodeMovementThresholdFraction
                                           * nodeMovementThresholdFraction ),
     bounceBeforeLastPath( functionValueForNanInput ),
     lastPathFalseSideNode( potentialFunction.NumberOfFieldVariables() ),
     lastPathTrueSideNode( potentialFunction.NumberOfFieldVariables() ),
     nodesConverged( false )
  {
    // This constructor is just an initialization list.
    // The factor of 0.25 in nodeMovementThresholdFractionSquared is to account
    // for nodeMovementThresholdFraction being given in reference to the length
    // of a single segment, while nodeMovementThresholdFractionSquared is used
    // in comparison with the square of the length of the difference of nodes
    // which are 2 segments apart.
  }

  MinuitOnPotentialPerpendicularToPath::~MinuitOnPotentialPerpendicularToPath()
  {
    // This does nothing.
  }


  // This minimizes the potential on hyperplanes which are as close to
  // perpendicular to the path given by lastPath, at numberOfVaryingNodes
  // equally-space points along that path, and the minima on those
  // hyperplanes are then used to construct the returned path. If the minima
  // are all sufficiently close to their starting points from lastPath,
  // nodesConverged is set to ture. (The bubble profile from the last path is
  // ignored, but there is an empty hook in the loop to allow derived classes
  // to use it.)
  TunnelPath const* MinuitOnPotentialPerpendicularToPath::TryToImprovePath(
                                                    TunnelPath const& lastPath,
                                      BubbleProfile const& bubbleFromLastPath )
  {
    // We assume that the nodes will converge with this call, and note if they
    // don't.
    lastPath.PutOnPathAt( currentHyperplaneOrigin,
                          segmentAuxiliaryLength );
    lastPath.PutOnPathAt( lastPathTrueSideNode,
                         ( segmentAuxiliaryLength + segmentAuxiliaryLength ) );

    // For the first varying node, we don't bother copying pathNodes.front()
    // into lastPathFalseSideNode, as it'll just be over-written by the first
    // iteration of the loop.
    SetParallelVector( returnPathNodes.front(),
                       lastPathTrueSideNode );
    SetCurrentMinuitSteps( 0.5 );
    // The starting step sizes for Minuit2 are set to be half the Euclidean
    // length of the difference between lastPathFalseSideNode and
    // lastPathTrueSideNode, times whatever internal factor Minuit2 uses (seems
    // to be 0.1 or 0.01), but the sizes are adapted as the minimization
    // proceeds anyway.
    SetUpHouseholderReflection();

    AccountForBubbleProfileAroundNode( 1,
                                       bubbleFromLastPath );
    RunMigradAndPutTransformedResultIn( returnPathNodes[ 1 ] );

    for( size_t nodeIndex( 2 );
         nodeIndex < numberOfVaryingNodes;
         ++nodeIndex )
    {
      lastPathFalseSideNode = currentHyperplaneOrigin;
      currentHyperplaneOrigin = lastPathTrueSideNode;
      lastPath.PutOnPathAt( lastPathTrueSideNode,
                            ( ( nodeIndex + 1 ) * segmentAuxiliaryLength ) );
      SetParallelVector( lastPathFalseSideNode,
                         lastPathTrueSideNode );
      double differenceSquared( 0.0 );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        differenceSquared += ( currentParallelComponent[ fieldIndex ]
                               * currentParallelComponent[ fieldIndex ] );
      }
      SetCurrentMinuitSteps( 0.5 );
      SetUpHouseholderReflection();
      AccountForBubbleProfileAroundNode( nodeIndex,
                                         bubbleFromLastPath );
      RunMigradAndPutTransformedResultIn( returnPathNodes[ nodeIndex ] );
      // This function leaves the last Minuit parameterization in
      // minuitResultAsUntransformedVector, so we can use it to determine "how
      // far the node rolled".
      if( minuitResultAsUntransformedVector.squaredNorm()
          > ( nodeMovementThresholdFractionSquared * differenceSquared ) )
      {
        nodesConverged = false;
      }
    }
    // Now we do the node just before the true vacuum. We don't bother copying
    // around all the nodes even though the code would be a bit easier to read.
    SetParallelVector( currentHyperplaneOrigin,
                       returnPathNodes.back() );
    currentHyperplaneOrigin = lastPathTrueSideNode;
    SetCurrentMinuitSteps( 0.5 );
    SetUpHouseholderReflection();
    AccountForBubbleProfileAroundNode( numberOfVaryingNodes,
                                       bubbleFromLastPath );
    RunMigradAndPutTransformedResultIn(
                                     returnPathNodes[ numberOfVaryingNodes ] );
    return new LinearSplineThroughNodes( returnPathNodes,
                                         nodeZeroParameterization,
                                         pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
