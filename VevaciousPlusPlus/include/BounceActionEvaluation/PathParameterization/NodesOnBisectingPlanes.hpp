/*
 * NodesOnBisectingPlanes.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef NODESONBISECTINGPLANES_HPP_
#define NODESONBISECTINGPLANES_HPP_

#include "CommonIncludes.hpp"
#include "Eigen/Dense"
#include "NodesOnPlanes.hpp"

namespace VevaciousPlusPlus
{

  class NodesOnBisectingPlanes : public NodesOnPlanes
  {
  public:
    NodesOnBisectingPlanes( size_t const numberOfFields,
                            size_t const numberOfIntermediateNodes );
    virtual ~NodesOnBisectingPlanes();


  protected:
    std::vector< size_t > adjustmentOrder;
    std::vector< std::pair< size_t, size_t > > sideNodeIndices;
    std::vector< Eigen::MatrixXd > rotationMatrices;


    // This just assumes that the nodes can be adjusted in the order which they
    // are visited in along the path, but can be over-ridden in derived
    // classes. (It also calls ValidAdjustmentIndex before returning the
    // index.)
    virtual size_t
    PathIndexFromAdjustmentIndex( size_t const adjustmentOrderIndex ) const;


    // This should add the perpendicular component from the parameterization
    // given by nodeParameterization along with startNode and endNode to
    // nodeVector.
    virtual void AddTransformedNode( std::vector< double >& nodeVector,
                                     std::vector< double > const& startNode,
                                     std::vector< double > const& endNode,
                     std::vector< double > const& nodeParameterization ) const;

    // This should return the false-vacuum-side node of the pair of nodes
    // from which the node at adjustmentOrderIndex should be set.
    virtual std::vector< double > const&
    FalseSideNode( size_t const adjustmentOrderIndex,
                   std::vector< std::vector< double > > const& nodeSet ) const
    { return nodeSet[ sideNodeIndices[ adjustmentOrderIndex ].first ]; }

    // This should return the true-vacuum-side node of the pair of nodes
    // from which the node at nodeIndex should be set.
    virtual std::vector< double > const&
    TrueSideNode( size_t const adjustmentOrderIndex,
                  std::vector< std::vector< double > > const& nodeSet ) const
    { return nodeSet[ sideNodeIndices[ adjustmentOrderIndex].second ]; }

    // This should return the fraction along the node difference vector that
    // the rotated plane should be shifted appropriate for
    // pathNodes[ nodeIndex ].
    virtual double ShiftFraction( size_t const nodeIndex ) const{ return 0.5; }
  };




  // This just assumes that the nodes can be adjusted in the order which they
  // are visited in along the path, but can be over-ridden in derived
  // classes.
  inline size_t NodesOnBisectingPlanes::PathIndexFromAdjustmentIndex(
                                      size_t const adjustmentOrderIndex ) const
  {
    ValidAdjustmentIndex( adjustmentOrderIndex );
    return adjustmentOrder[ adjustmentOrderIndex ];
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESONBISECTINGPLANES_HPP_ */
