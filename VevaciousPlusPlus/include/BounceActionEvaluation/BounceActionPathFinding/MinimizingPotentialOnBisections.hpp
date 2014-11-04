/*
 * MinimizingPotentialOnBisections.hpp
 *
 *  Created on: Sep 11, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINIMIZINGPOTENTIALONBISECTIONS_HPP_
#define MINIMIZINGPOTENTIALONBISECTIONS_HPP_

#include "CommonIncludes.hpp"
#include "Eigen/Dense"
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
                                     unsigned int const minuitStrategy = 1,
                                    double const minuitToleranceFraction = 0.5,
                                     size_t const maximumNumberOfNodes = 254 );
    virtual ~MinimizingPotentialOnBisections();

    needs number of uphill improvements!

    need to incorporate use of currentParallelComponent properly!

    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node.
    virtual double
    operator()( std::vector< double > const& nodeParameterization ) const;


  protected:
    size_t const maximumNumberOfNodes;
    bool inSegmentSplittingStage;
    size_t finalNumberOfNodes;
    double segmentAuxiliaryLength;
    std::vector< double > bisectionOrigin;
    std::vector< std::vector< double > > lastNodes;


    // This should internally change the derived class's strategy based on
    // whether or not the last attempt to improve the path lowered the bounce
    // action successfully or not, and the return true if there are still
    // improvements to try, false if there is no improvement to try.
    virtual bool NodesCanStillBeImproved( bool const lastMoveDidImprove );

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
    virtual void TryToImproveNodes();

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

    // This inserts a new node in between each pair of old nodes in pathNodes
    // where the potential in the bisecting hyperplane is minimized.
    bool SplitSegments();

    // This forms a new path from splines through pathNodes, then adjusts the
    // nodes of pathNodes based on nodes chosen equally spaced along the path.
    void AdjustNodes();
  };




  // It may seem unwise to have this object call Minuit on itself, but really
  // it's just a handy way of keeping the minimization function within the
  // class that ends up finding its minimum. In this case,
  // nodeParameterization is just the parameterization of the node.
  inline double MinimizingPotentialOnBisections::operator()(
                      std::vector< double > const& nodeParameterization ) const
  {
    Eigen::VectorXd const transformedNode( reflectionMatrix
                                 * UntransformedNode( nodeParameterization ) );
    std::vector< double > fieldConfiguration( bisectionOrigin );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      fieldConfiguration[ fieldIndex ] += transformedNode( fieldIndex );
    }
    return potentialFunction( fieldConfiguration,
                              pathTemperature );
  }

  // If this MinimizingPotentialOnBisections is still in the segment-splitting
  // stage, it either returns true if the last move did lower the bounce action
  // or it reverts to the last set of nodes, noting that it is no longer in the
  // segment-splitting stage, and then returns true. If it is not in the
  // segment-splitting stage when this function is called, then it returns true
  // if the last move did lower the bounce action or false if it didn't, as
  // iterating the improvement after it gets worse the first time seems to just
  // usually make things worse without ever recovering.
  inline bool MinimizingPotentialOnBisections::NodesCanStillBeImproved(
                                                bool const lastMoveDidImprove )
  {
    if( inSegmentSplittingStage )
    {
      if( !lastMoveDidImprove
          &&
          ( pathNodes.size() > 3 ) )
      {
        pathNodes = lastNodes;
        inSegmentSplittingStage = false;
        finalNumberOfNodes = pathNodes.size();
        segmentAuxiliaryLength
        = ( 1.0 / static_cast< double >( finalNumberOfNodes - 1 ) );
      }
      return true;
    }
    else
    {
      return lastMoveDidImprove;
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
    lastNodes = pathNodes;
    inSegmentSplittingStage = true;
  }

  // This inserts a new node in between each pair of old nodes in pathNodes
  // where the potential in the bisecting hyperplane is minimized. It returns
  // true if the number of nodes between the vacua is now at least
  // minimumNumberOfNodes.
  inline bool MinimizingPotentialOnBisections::SplitSegments()
  {
    lastNodes = pathNodes;
    pathNodes.resize( ( ( 2 * lastNodes.size() ) - 1 ),
                      lastNodes.back() );
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
    // We return whether another split would still give less nodes than the
    // maximum allowed.
    return ( ( 2 * ( pathNodes.size() - 2 ) ) < maximumNumberOfNodes );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGPOTENTIALONBISECTIONS_HPP_ */
