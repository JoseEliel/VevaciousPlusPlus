/*
 * MinimizingPotentialOnBisections.cpp
 *
 *  Created on: Sep 11, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinimizingPotentialOnBisections.hpp"

namespace VevaciousPlusPlus
{

  MinimizingPotentialOnBisections::MinimizingPotentialOnBisections(
                                    PotentialFunction const& potentialFunction,
                                               MinuitBetweenPaths* pathRefiner,
                                             size_t const minimumNumberOfNodes,
                                            double const moveToleranceFraction,
                                              size_t const movesPerImprovement,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinimizingPotentialOnHypersurfaces( potentialFunction,
                                        pathRefiner,
                                        movesPerImprovement,
                                        minuitStrategy,
                                        minuitToleranceFraction ),
    minimumNumberOfNodes( minimumNumberOfNodes ),
    moveToleranceFractionSquared( moveToleranceFraction
                                  * moveToleranceFraction ),
    inSegmentSplittingStage( true ),
    finalNumberOfNodes( 0 ),
    segmentAuxiliaryLength( NAN )
  {
    // This constructor is just an initialization list.
  }

  MinimizingPotentialOnBisections::~MinimizingPotentialOnBisections()
  {
    // This does nothing.
  }


  // This sets up currentParallelComponent to be the vector from
  // falseSideNode to trueSideNode, and sets reflectionMatrix appropriately,
  // and sets bisectionOrigin to be halfway between falseSideNode and
  // trueSideNode, then sets nodeToSet to be the point on the bisecting
  // hyperplane where the potential is minimized.
  void MinimizingPotentialOnBisections::PlaceOnMinimumInBisection(
                                              std::vector< double >& nodeToSet,
                                    std::vector< double > const& falseSideNode,
                                    std::vector< double > const& trueSideNode )
  {
    SetParallelVector( falseSideNode,
                       trueSideNode );
    SetUpHouseholderReflection();
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      bisectionOrigin[ fieldIndex ] = ( 0.5 * ( falseSideNode[ fieldIndex ]
                                              + trueSideNode[ fieldIndex ] ) );
    }
    ROOT::Minuit2::MnMigrad mnMigrad( *this,
                                      nodeZeroParameterization,
                                      minuitInitialSteps,
                                      minuitStrategy );
    MinuitMinimum minuitResult( ( numberOfFields - 1 ),
                                mnMigrad( 0,
                                          currentMinuitTolerance ) );
    Eigen::VectorXd const transformedNode( reflectionMatrix
                        * UntransformedNode( minuitResult.VariableValues() ) );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      nodeToSet[ fieldIndex ] = ( bisectionOrigin[ fieldIndex ]
                                  + transformedNode( fieldIndex ) );
    }
  }

  // This forms a new path from splines through pathNodes, then adjusts the
  // nodes of pathNodes based on nodes chosen equally spaced along the path.
  void MinimizingPotentialOnBisections::AdjustNodes()
  {
    nodesConverged = true;
    LinearSplineThroughNodes lastPath( pathNodes,
                                       std::vector< double >(),
                                       pathTemperature );
    std::vector< double > nextNode( numberOfFields );
    std::vector< double > oldPosition( numberOfFields );
    for( size_t nodeIndex( 1 );
         nodeIndex < ( finalNumberOfNodes - 2 );
         ++nodeIndex )
    {
      lastPath.PutOnPathAt( oldPosition,
                            ( nodeIndex * segmentAuxiliaryLength ) );
      lastPath.PutOnPathAt( nextNode,
                            ( ( nodeIndex + 1 ) * segmentAuxiliaryLength ) );
      PlaceOnMinimumInBisection( pathNodes[ nodeIndex ],
                                 pathNodes[ nodeIndex - 1 ],
                                 nextNode );
      if( NodeMovedBeyondThreshold( oldPosition,
                                    pathNodes[ nodeIndex ],
                                    pathNodes[ nodeIndex - 1 ] ) )
      {
        nodesConverged = false;
      }
    }
    lastPath.PutOnPathAt( oldPosition,
                          ( 1.0 - segmentAuxiliaryLength ) );
    PlaceOnMinimumInBisection( pathNodes[ finalNumberOfNodes - 2 ],
                               pathNodes[ finalNumberOfNodes - 3 ],
                               pathNodes.back() );
    if( NodeMovedBeyondThreshold( oldPosition,
                                  pathNodes[ finalNumberOfNodes - 2 ],
                                  pathNodes[ finalNumberOfNodes - 3 ] ) )
    {
      nodesConverged = false;
    }
  }

} /* namespace VevaciousPlusPlus */
