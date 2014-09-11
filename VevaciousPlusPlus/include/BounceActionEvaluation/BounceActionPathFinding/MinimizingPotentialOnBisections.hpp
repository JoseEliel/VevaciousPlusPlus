/*
 * MinimizingPotentialOnBisections.hpp
 *
 *  Created on: Sep 11, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINIMIZINGPOTENTIALONBISECTIONS_HPP_
#define MINIMIZINGPOTENTIALONBISECTIONS_HPP_

#include "CommonIncludes.hpp"
#include "MinimizingPotentialOnHypersurfaces.hpp"
#include "../PathParameterization/LinearSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class MinimizingPotentialOnBisections :
                                      public MinimizingPotentialOnHypersurfaces
  {
  public:
    MinimizingPotentialOnBisections(
                                    PotentialFunction const& potentialFunction,
                                     MinuitBetweenPaths* pathRefiner,
                                     size_t const minimumNumberOfNodes,
                                     double const moveToleranceFraction,
                                     size_t const movesPerImprovement = 100,
                                     unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5 );
    virtual ~MinimizingPotentialOnBisections();


    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node.
    virtual double
    operator()( std::vector< double > const& nodeParameterization ) const;


  protected:
    size_t const minimumNumberOfNodes;
    double const moveToleranceFractionSquared;
    bool inSegmentSplittingStage;
    size_t finalNumberOfNodes;
    double segmentAuxiliaryLength;
    std::vector< double > bisectionOrigin;


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
    virtual void ImproveNodes();

    // This sets pathNodes to just be the false vacuum node and the true vacuum
    // node, and also sets nodesConverged and inSegmentSplittingStage
    // appropriately.
    virtual void SetNodesForInitialPath( PotentialMinimum const& falseVacuum,
                                         PotentialMinimum const& trueVacuum );

    // This sets up currentParallelComponent to be the vector from
    // falseSideNode to trueSideNode, and sets reflectionMatrix appropriately,
    // and sets bisectionOrigin to be halfway between falseSideNode and
    // trueSideNode, then sets nodeToSet to be the point on the bisecting
    // hyperplane where the potential is minimized.
    void PlaceOnMinimumInBisection( std::vector< double >& nodeToSet,
                                    std::vector< double > const& falseSideNode,
                                   std::vector< double > const& trueSideNode );

    // This returns true if the node moved more than
    // ( moveToleranceFractionSquared )^(1/2) times the distance between where
    // it originally was and its false-side nearest neighbor.
    bool NodeMovedBeyondThreshold( std::vector< double > const& newPosition,
                                   std::vector< double > const& oldPosition,
                            std::vector< double > const& falseSideNode ) const;

    // This inserts a new node in between each pair of old nodes in pathNodes
    // where the potential in the bisecting hyperplane is minimized.
    bool SplitSegments();

    // This forms a new path from splines through pathNodes, then adjusts the
    // nodes of pathNodes based on nodes chosen equally spaced along the path.
    void AdjustNodes();
  };


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
  inline void MinimizingPotentialOnBisections::ImproveNodes()
  {
    if( inSegmentSplittingStage )
    {
      inSegmentSplittingStage = SplitSegments();
      if( !inSegmentSplittingStage )
      {
        finalNumberOfNodes = pathNodes.size();
        segmentAuxiliaryLength
        = ( 1.0 / static_cast< double >( finalNumberOfNodes - 1 ) );
      }
    }
    else
    {
      AdjustNodes();
    }
  }

  // This sets pathNodes to just be the false vacuum node and the true vacuum
  // node, and also sets nodesConverged and inSegmentSplittingStage
  // appropriately.
  inline void MinimizingPotentialOnBisections::SetNodesForInitialPath(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    pathNodes.resize( 2 );
    pathNodes.front() = falseVacuum.FieldConfiguration();
    pathNodes.back() = trueVacuum.FieldConfiguration();
    nodesConverged = false;
    inSegmentSplittingStage = true;
  }

  // This returns true if the node moved more than
  // ( moveToleranceFractionSquared )^(1/2) times the distance between where it
  // originally was and its false-side nearest neighbor.
  inline bool MinimizingPotentialOnBisections::NodeMovedBeyondThreshold(
                                      std::vector< double > const& newPosition,
                                      std::vector< double > const& oldPosition,
                             std::vector< double > const& falseSideNode ) const
  {
    double displacementSquared( 0.0 );
    double comparisonSquared( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double
      fieldDifference( newPosition[ fieldIndex ] - oldPosition[ fieldIndex ] );
      displacementSquared += ( fieldDifference * fieldDifference );
      fieldDifference
      = ( oldPosition[ fieldIndex ] - falseSideNode[ fieldIndex ] );
      comparisonSquared += ( fieldDifference * fieldDifference );
    }
    return ( displacementSquared > ( comparisonSquared
                                     * moveToleranceFractionSquared ) );
  }

  // This inserts a new node in between each pair of old nodes in pathNodes
  // where the potential in the bisecting hyperplane is minimized. It returns
  // true if the number of nodes between the vacua is now at least
  // minimumNumberOfNodes.
  inline bool MinimizingPotentialOnBisections::SplitSegments()
  {
    std::vector< std::vector< double > > const lastNodes( pathNodes );
    pathNodes.resize( ( 2 * lastNodes.size() ) - 1 );
    // pathNodes.front() = lastNodes.front();
    // This line is redundant because std::vector::resize( ... ) should not
    // change what is in pathNodes up to the new size.
    for( size_t nodeIndex( 1 );
         nodeIndex < lastNodes.size();
         ++nodeIndex )
    {
      pathNodes[ 2 * nodeIndex ] = lastNodes[ nodeIndex ];
      PlaceOnMinimumInBisection( pathNodes[ ( 2 * nodeIndex ) - 1 ],
                                 lastNodes[ nodeIndex - 1 ],
                                 lastNodes[ nodeIndex ] );
    }
    return ( ( pathNodes.size() - 2 ) >= minimumNumberOfNodes );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGPOTENTIALONBISECTIONS_HPP_ */
