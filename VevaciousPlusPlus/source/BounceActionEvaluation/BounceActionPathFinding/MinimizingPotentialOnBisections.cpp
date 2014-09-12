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
    segmentAuxiliaryLength( NAN ),
    bisectionOrigin( numberOfFields )
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
    SetCurrentMinuitSteps();
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      bisectionOrigin[ fieldIndex ] = ( 0.5 * ( falseSideNode[ fieldIndex ]
                                              + trueSideNode[ fieldIndex ] ) );
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
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
    std::cout << std::endl;/**/

    ROOT::Minuit2::MnMigrad mnMigrad( *this,
                                      nodeZeroParameterization,
                                      minuitInitialSteps,
                                      minuitStrategy );
    MinuitMinimum minuitResult( ( numberOfFields - 1 ),
                                mnMigrad( 0,
                                          currentMinuitTolerance ) );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "minuitResult =" << std::endl << minuitResult.AsDebuggingString();
    std::cout << std::endl;/**/

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
    /**/std::cout << std::endl << "debugging:"
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
    std::cout << std::endl;/**/
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

  // This inserts a new node in between each pair of old nodes in pathNodes
  // where the potential in the bisecting hyperplane is minimized. It returns
  // true if the number of nodes between the vacua is now at least
  // minimumNumberOfNodes.
  bool MinimizingPotentialOnBisections::SplitSegments()
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "MinimizingPotentialOnBisections::SplitSegments() called." << std::endl;
    std::cout << "pathNodes = { { ";
    for( std::vector< std::vector< double > >::const_iterator
         pathNode( pathNodes.begin() );
         pathNode < pathNodes.end();
         ++pathNode )
    {
      if( pathNode > pathNodes.begin() )
      {
        std::cout << " }," << std::endl << "{ ";
      }
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << (*pathNode)[ fieldIndex ];
      }
    }
    std::cout << " } }";
    std::cout << std::endl;/**/

    std::vector< std::vector< double > > const lastNodes( pathNodes );
    pathNodes.resize( ( ( 2 * lastNodes.size() ) - 1 ),
                      lastNodes.back() );
    // pathNodes.front() = lastNodes.front();
    // This line is redundant because std::vector::resize( ... ) should not
    // change what is in pathNodes up to the new size.
    for( size_t nodeIndex( 1 );
         nodeIndex < lastNodes.size();
         ++nodeIndex )
    {
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "nodeIndex = " << nodeIndex;
      std::cout << " pathNodes = { { ";
      for( std::vector< std::vector< double > >::const_iterator
           pathNode( pathNodes.begin() );
           pathNode < pathNodes.end();
           ++pathNode )
      {
        if( pathNode > pathNodes.begin() )
        {
          std::cout << " }," << std::endl << "{ ";
        }
        for( size_t fieldIndex( 0 );
             fieldIndex < pathNode->size();
             ++fieldIndex )
        {
          if( fieldIndex > 0 )
          {
            std::cout << ", ";
          }
          std::cout << (*pathNode)[ fieldIndex ];
        }
      }
      std::cout << " } }";
      std::cout << std::endl;/**/

      pathNodes[ 2 * nodeIndex ] = lastNodes[ nodeIndex ];
      PlaceOnMinimumInBisection( pathNodes[ ( 2 * nodeIndex ) - 1 ],
                                 lastNodes[ nodeIndex - 1 ],
                                 lastNodes[ nodeIndex ] );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl << "now pathNodes = { { ";
      for( std::vector< std::vector< double > >::const_iterator
           pathNode( pathNodes.begin() );
           pathNode < pathNodes.end();
           ++pathNode )
      {
        if( pathNode > pathNodes.begin() )
        {
          std::cout << " }," << std::endl << "{ ";
        }
        for( size_t fieldIndex( 0 );
             fieldIndex < pathNode->size();
             ++fieldIndex )
        {
          if( fieldIndex > 0 )
          {
            std::cout << ", ";
          }
          std::cout << (*pathNode)[ fieldIndex ];
        }
      }
      std::cout << " } }";
      std::cout << std::endl;/**/
    }
    return ( ( pathNodes.size() - 2 ) < minimumNumberOfNodes );
  }

} /* namespace VevaciousPlusPlus */
