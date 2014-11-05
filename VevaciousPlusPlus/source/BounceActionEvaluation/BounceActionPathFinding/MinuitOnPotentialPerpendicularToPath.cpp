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
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
     MinuitOnHypersurfaces( potentialFunction,
                            numberOfPathSegments,
                            minuitStrategy,
                            minuitToleranceFraction ),
     numberOfAllowedWorsenings( numberOfAllowedWorsenings ),
     bounceBeforeLastPath( functionValueForNanInput ),
     lastPathFalseSideNode( potentialFunction.NumberOfFieldVariables() ),
     lastPathTrueSideNode( potentialFunction.NumberOfFieldVariables() )
  {
    // This constructor is just an initialization list.
  }

  MinuitOnPotentialPerpendicularToPath::~MinuitOnPotentialPerpendicularToPath()
  {
    // This does nothing.
  }


  // This takes the bounce action from bubbleFromLastPath and updates
  // numberOfAllowedWorsenings based on whether it was an improvement on the
  // previous path. Then it returns false if too many paths have been tried
  // which just made the action bigger, true otherwise.
  bool MinuitOnPotentialPerpendicularToPath::PathCanBeImproved(
                                      BubbleProfile const& bubbleFromLastPath )
  {
    if( bubbleFromLastPath.BounceAction() >= bounceBeforeLastPath )
    {
      --numberOfAllowedWorsenings;
    }
    bounceBeforeLastPath = bubbleFromLastPath.BounceAction();
    return ( numberOfAllowedWorsenings > 0 );
  }

  // This minimizes the potential on hyperplanes which are as close to
  // perpendicular to the path given by lastPath, at numberOfVaryingNodes
  // equally-space points along that path, and the minima on those
  // hyperplanes are then used to construct the returned path. (The bubble
  // profile from the last path is ignored, but there is an empty hook in the
  // loop to allow derived classes to use it.)
  TunnelPath const* MinuitOnPotentialPerpendicularToPath::TryToImprovePath(
                                                    TunnelPath const& lastPath,
                                      BubbleProfile const& bubbleFromLastPath )
  {
    lastPath.PutOnPathAt( currentHyperplaneOrigin,
                          segmentAuxiliaryLength );
    lastPath.PutOnPathAt( lastPathTrueSideNode,
                         ( segmentAuxiliaryLength + segmentAuxiliaryLength ) );

    // For the first varying node, we don't bother copying pathNodes.front()
    // into lastPathFalseSideNode, as it'll just be over-written by the first
    // iteration of the loop.
    SetParallelVector( pathNodes.front(),
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
    RunMigradAndPutTransformedResultIn( pathNodes[ 1 ] );

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
      SetCurrentMinuitSteps( 0.5 );
      SetUpHouseholderReflection();
      AccountForBubbleProfileAroundNode( nodeIndex,
                                         bubbleFromLastPath );
      RunMigradAndPutTransformedResultIn( pathNodes[ nodeIndex ] );
    }
    // Now we do the node just before the true vacuum. We don't bother copying
    // around all the nodes even though the code would be a bit easier to read.
    SetParallelVector( currentHyperplaneOrigin,
                       pathNodes.back() );
    currentHyperplaneOrigin = lastPathTrueSideNode;
    SetCurrentMinuitSteps( 0.5 );
    SetUpHouseholderReflection();
    AccountForBubbleProfileAroundNode( numberOfVaryingNodes,
                                       bubbleFromLastPath );
    RunMigradAndPutTransformedResultIn( pathNodes[ numberOfVaryingNodes ] );
    return new LinearSplineThroughNodes( pathNodes,
                                         nodeZeroParameterization,
                                         pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
