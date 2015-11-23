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
                                             size_t const numberOfPathSegments,
                                           int const numberOfAllowedWorsenings,
                                    double const nodeMovementThresholdFraction,
                                                    double const dampingFactor,
                       std::vector< double > const neighborDisplacementWeights,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
     MinuitOnHypersurfaces( numberOfPathSegments,
                            minuitStrategy,
                            minuitToleranceFraction ),
     numberOfAllowedWorsenings( numberOfAllowedWorsenings ),
     numberOfWorseningsSoFar( 0 ),
     nodeMovementThresholdFractionSquared( 0.25 * nodeMovementThresholdFraction
                                           * nodeMovementThresholdFraction ),
     dampingFactor( dampingFactor ),
     neighborDisplacementWeights( neighborDisplacementWeights ),
     lastPathNodes( ( numberOfPathSegments + 1 ),
                    Eigen::VectorXd::Zero( numberOfFields ) ),
     nodeDisplacements( lastPathNodes ),
     bounceBeforeLastPath( functionValueForNanInput ),
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
  // nodesConverged is set to true. (The bubble profile from the last path is
  // ignored, but there is an empty hook in the loop to allow derived classes
  // to use it.)
  TunnelPath const* MinuitOnPotentialPerpendicularToPath::TryToImprovePath(
                                                    TunnelPath const& lastPath,
                                      BubbleProfile const& bubbleFromLastPath )
  {
    // First we set up the nodes corresponding to the last path, ignoring the
    // ends, which should have been set correctly by SetNodesForInitialPath.
    for( size_t nodeIndex( 1 );
         nodeIndex <= numberOfVaryingNodes;
         ++nodeIndex )
    {
      lastPath.PutOnPathAt( lastPathNodes[ nodeIndex ],
                            ( nodeIndex * segmentAuxiliaryLength ) );
    }

    // We assume that the nodes will converge with this call, and note if they
    // don't.
    nodesConverged = true;

    // Next we calculate and store the displacements, ignoring the ends, which
    // should not get displaced.
    for( size_t nodeIndex( 1 );
         nodeIndex <= numberOfVaryingNodes;
         ++nodeIndex )
    {
      currentHyperplaneOrigin = lastPathNodes[ nodeIndex ];
      SetParallelVector( lastPathNodes[ nodeIndex - 1 ],
                         lastPathNodes[ nodeIndex + 1 ] );
      SetCurrentMinuitSteps( 0.1 );
      // The starting step sizes for Minuit2 are set to be a tenth of the
      // Euclidean length of the difference between
      // lastPathNodes[ nodeIndex - 1 ] and lastPathNodes[ nodeIndex + 1 ],
      // times whatever internal factor Minuit2 uses (seems to be 0.1 or 0.01),
      // but the sizes are adapted as the minimization proceeds anyway.
      SetUpHouseholderReflection();
      AccountForBubbleProfileAroundNode( nodeIndex,
                                         bubbleFromLastPath );
      nodeDisplacements[ nodeIndex ] = RunMigradAndReturnDisplacement();
      if( nodeDisplacements[ nodeIndex ].squaredNorm()
          > ( nodeMovementThresholdFractionSquared
              * currentParallelComponent.squaredNorm() ) )
      {
        nodesConverged = false;
      }
    }

    // Now we set returnPathNodes based on lastPathNodes plus scaled weighted
    // averages of nodeDisplacements, ignoring the ends, which should not get
    // displaced.
    for( size_t nodeIndex( 1 );
         nodeIndex <= numberOfVaryingNodes;
         ++nodeIndex )
    {
      // It might be more efficient to average individual components and put
      // them directly into returnPathNodes[ nodeIndex ], but this way is
      // clearer. It might get optimized to be as efficient as doing individual
      // components anyway.
      Eigen::VectorXd
      weightedAverageDisplacement( nodeDisplacements[ nodeIndex ] );
      double weightNormalization( 1.0 );
      size_t const weightsSize( neighborDisplacementWeights.size() );
      // We use neighborIndex going from 1 to weightsSize because it's]
      // conceptually easier to think of number of steps to get to the
      // neighbor. Hence the nearest neighbor has neighborIndex of 1, and
      // the nearest neighbor displacements for node nodeIndex are
      // nodeDisplacements[ nodeIndex - 1 ] and
      // nodeDisplacements[ nodeIndex + 1 ], making things much easier to think
      // about, bearing in mind that the weight is found in
      // neighborDisplacementWeights[ neighborIndex - 1 ].
      for( size_t neighborIndex( 1 );
           neighborIndex <= weightsSize;
           ++neighborIndex )
      {
        double const
        neighborWeight( neighborDisplacementWeights[ neighborIndex - 1 ] );
        if( nodeIndex > neighborIndex )
        {
          // We average in the false-side neighbors if they are not at the end
          // or beyond.
          weightNormalization += neighborWeight;
          weightedAverageDisplacement += ( neighborWeight
                            * nodeDisplacements[ nodeIndex - neighborIndex ] );
        }
        if( ( nodeIndex + neighborIndex ) <= numberOfVaryingNodes )
        {
          // We average in the true-side neighbors if they are not at the end
          // or beyond.
          weightNormalization += neighborWeight;
          weightedAverageDisplacement += ( neighborWeight
                            * nodeDisplacements[ nodeIndex + neighborIndex ] );
        }
      }

      double const displacementScaling( dampingFactor / weightNormalization );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        returnPathNodes[ nodeIndex ][ fieldIndex ]
        = ( lastPathNodes[ nodeIndex ]( fieldIndex )
            + ( displacementScaling
                * weightedAverageDisplacement( fieldIndex ) ) );
      }
    }

    // Finally we make a new path through returnPathNodes.
    return new LinearSplineThroughNodes( returnPathNodes,
                                         nodeZeroParameterization,
                                         pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
