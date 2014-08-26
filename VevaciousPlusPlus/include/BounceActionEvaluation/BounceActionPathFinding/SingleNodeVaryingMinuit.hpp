/*
 * SingleNodeVaryingMinuit.hpp
 *
 *  Created on: Jul 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SINGLENODEVARYINGMINUIT_HPP_
#define SINGLENODEVARYINGMINUIT_HPP_

#include "CommonIncludes.hpp"
#include "MinuitPathFinder.hpp"
#include "Minuit2/MnMigrad.h"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitMinimum.hpp"
#include "../PathParameterization/PathFromNodesFactory.hpp"
#include "../PathParameterization/NodesFromParameterization.hpp"

namespace VevaciousPlusPlus
{

  class SingleNodeVaryingMinuit : public MinuitPathFinder
  {
  public:
    SingleNodeVaryingMinuit( PotentialFunction const& potentialFunction,
                             PathFromNodesFactory* const pathFactory,
                             size_t const movesPerImprovement = 100,
                             unsigned int const minuitStrategy = 1,
                             double const minuitToleranceFraction = 0.5,
                             double const nodeMoveThreshold = 0.01 );
    virtual ~SingleNodeVaryingMinuit();


    // This resets the BouncePathFinder so that it sets up currentPath as its
    // initial path between the given vacua. It also resets pathCanBeImproved
    // and sets pathTemperature appropriately.
    virtual void SetInitialPath( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 TunnelPath const* startingPath = NULL,
                                 double const pathTemperature = 0.0 );

    // This allows Minuit2 to move each node a set number of times to try to
    // minimize the potential at that node, and then sets the path from the set
    // of nodes.
    virtual void ImprovePath();


  protected:
    PotentialFunction const& potentialFunction;
    PathFromNodesFactory* const pathFactory;
    NodesFromParameterization& pathNodes;
    std::vector< MinuitMinimum > currentMinuitResults;
    size_t currentNodeIndex;
    double const nodeMoveThresholdSquared;


    // This returns true if the distance between minuitNode and the node at
    // currentNodeIndex (before updating the node) is greater than a threshold
    // fraction times the smaller of the Euclidean lengths between the node at
    // currentNodeIndex and its 2 nearest neighbors.
    bool NodeMovedMoreThanThreshold(
                               std::vector< double > const& minuitNode ) const;

    // This returns the sum of the squares of the differences of the elements
    // of firstVector and secondVector.
    double
    DifferenceEuclideanLengthSquared( std::vector< double > const& firstVector,
                             std::vector< double > const& secondVector ) const;
  };




  // This returns true if the distance between minuitNode and the node at
  // currentNodeIndex (before updating the node) is greater than a threshold
  // fraction times the smaller of the Euclidean lengths between the node at
  // currentNodeIndex and its 2 nearest neighbors.
  inline bool SingleNodeVaryingMinuit::NodeMovedMoreThanThreshold(
                                std::vector< double > const& minuitNode ) const
  {
    std::vector< std::vector< double > > const&
    nodeVectors( pathNodes.PathNodes() );
    std::vector< double > const&
    currentNode( nodeVectors[ currentNodeIndex ] );
    double const nodeMoveSquared( DifferenceEuclideanLengthSquared( minuitNode,
                                                               currentNode ) );
    double const
    falseSideSquared( DifferenceEuclideanLengthSquared( currentNode,
                                       nodeVectors[ currentNodeIndex - 1 ] ) );
    double const
    trueSideSquared( DifferenceEuclideanLengthSquared( currentNode,
                                       nodeVectors[ currentNodeIndex + 1 ] ) );
    double thresholdSquared( nodeMoveThresholdSquared * falseSideSquared );
    if( trueSideSquared < falseSideSquared )
    {
      thresholdSquared = ( nodeMoveThresholdSquared * trueSideSquared );
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "SingleNodeVaryingMinuit::NodeMovedMoreThanThreshold( minuitNode = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < minuitNode.size();
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << minuitNode[ fieldIndex ];
    }
    std::cout << " } ) called. currentNodeIndex = " << currentNodeIndex
    << ", currentNode = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < currentNode.size();
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << currentNode[ fieldIndex ];
    }
    std::cout << " }, nodeVectors[ currentNodeIndex - 1 ] = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < nodeVectors[ currentNodeIndex - 1 ].size();
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << nodeVectors[ currentNodeIndex - 1 ][ fieldIndex ];
    }
    std::cout << " }, nodeVectors[ currentNodeIndex + 1 ] = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < nodeVectors[ currentNodeIndex + 1 ].size();
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << nodeVectors[ currentNodeIndex + 1 ][ fieldIndex ];
    }
    std::cout << " }, nodeMoveThresholdSquared = " << nodeMoveThresholdSquared
    << ", nodeMoveSquared = " << nodeMoveSquared << ", falseSideSquared = "
    << falseSideSquared << ", trueSideSquared = " << trueSideSquared
    << ", thresholdSquared = " << thresholdSquared
    << ", ( nodeMoveSquared > thresholdSquared ) = "
    << ( nodeMoveSquared > thresholdSquared );
    std::cout << std::endl;/**/

    return ( nodeMoveSquared > thresholdSquared );
  }

  // This returns the sum of the squares of the differences of the elements
  // of firstVector and secondVector.
  inline double SingleNodeVaryingMinuit::DifferenceEuclideanLengthSquared(
                                      std::vector< double > const& firstVector,
                              std::vector< double > const& secondVector ) const
  {
    double lengthSquared( 0.0 );
    for( size_t whichIndex( 0 );
         whichIndex < firstVector.size();
         ++whichIndex )
    {
      double const elementDifference( firstVector[ whichIndex ]
                                      - secondVector[ whichIndex ] );
      lengthSquared += ( elementDifference * elementDifference );
    }
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "SingleNodeVaryingMinuit::DifferenceEuclideanLengthSquared( firstVector"
    << " = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < firstVector.size();
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << firstVector[ fieldIndex ];
    }
    std::cout << " }, secondVector = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < secondVector.size();
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << secondVector[ fieldIndex ];
    }
    std::cout << " } ) about to return lengthSquared = " << lengthSquared;
    std::cout << std::endl;/**/
    return lengthSquared;
  }

} /* namespace VevaciousPlusPlus */
#endif /* SINGLENODEVARYINGMINUIT_HPP_ */
