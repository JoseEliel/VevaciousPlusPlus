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
                                             unsigned int const minuitStrategy,
                                          double const minuitToleranceFraction,
                                          size_t const maximumNumberOfNodes ) :
    MinimizingPotentialOnHypersurfaces( potentialFunction,
                                        minuitStrategy,
                                        minuitToleranceFraction ),
    maximumNumberOfNodes( maximumNumberOfNodes ),
    inSegmentSplittingStage( true ),
    finalNumberOfNodes( 0 ),
    segmentAuxiliaryLength( NAN ),
    bisectionOrigin( numberOfFields ),
    lastNodes()
  {
    // This constructor is just an initialization list.
  }

  MinimizingPotentialOnBisections::~MinimizingPotentialOnBisections()
  {
    // This does nothing.
  }


  // This sets up the nodes by minimizing the potential on hyperplanes
  // bisecting the vectors between pairs of nodes. In the first stage, the
  // segment between the vacuum is bisected, forming 2 segments, between the
  // false vacuum and the new node, and between the new node and the true
  // vacuum. Next those segments are themselves bisected, then the next set
  // of segments is bisected, and so on, until there are at least
  // minimumNumberOfNodes nodes between the vacua. When there are enough
  // nodes, the next stage begins, where the nodes form a path from splines,
  // then the nodes are equally spaced along this path, and each is moved in
  // the hyperplane bisecting the line between the node's nearest neighbors
  // to minimize the potential at that node in the hyperplane. This continues
  // until all of the nodes moved less than moveToleranceFraction times the
  // average segment length.
  void MinimizingPotentialOnBisections::TryToImproveNodes()
  {
    if( inSegmentSplittingStage )
    {
      inSegmentSplittingStage = SplitSegments();
    }
    else
    {
      AdjustNodes();
    }
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
    SetCurrentMinuitSteps();
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      bisectionOrigin[ fieldIndex ] = ( 0.5 * ( falseSideNode[ fieldIndex ]
                                              + trueSideNode[ fieldIndex ] ) );
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "MinimizingPotentialOnBisections::PlaceOnMinimumInBisection("
    << " nodeToSet = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << nodeToSet[ fieldIndex ];
    }
    std::cout << " }," << std::endl << "falseSideNode = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << falseSideNode[ fieldIndex ];
    }
    std::cout << " }," << std::endl << "trueSideNode = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << trueSideNode[ fieldIndex ];
    }
    std::cout << " } ) about to run Minuit2." << std::endl;
    std::cout << "bisectionOrigin = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << bisectionOrigin[ fieldIndex ];
    }
    std::cout << " }," << std::endl << "currentParallelComponent = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << currentParallelComponent[ fieldIndex ];
    }
    std::cout << " }";
    std::cout << std::endl;*/

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

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PlaceOnMinimumInBisection( ... ), finished, nodeToSet = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << nodeToSet[ fieldIndex ];
    }
    std::cout << " }";
    std::cout << std::endl;*/
  }

  // This forms a new path from splines through pathNodes, then adjusts the
  // nodes of pathNodes based on nodes chosen equally spaced along the path.
  void MinimizingPotentialOnBisections::AdjustNodes()
  {
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
    }
    lastPath.PutOnPathAt( oldPosition,
                          ( 1.0 - segmentAuxiliaryLength ) );
    PlaceOnMinimumInBisection( pathNodes[ finalNumberOfNodes - 2 ],
                               pathNodes[ finalNumberOfNodes - 3 ],
                               pathNodes.back() );
  }

} /* namespace VevaciousPlusPlus */
