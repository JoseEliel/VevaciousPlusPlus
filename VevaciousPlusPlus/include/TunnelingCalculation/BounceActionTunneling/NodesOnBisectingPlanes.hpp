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
    NodesOnBisectingPlanes( std::vector< double > const& falseVacuum,
                            std::vector< double > const& trueVacuum,
                            size_t const numberOfIntermediateNodes );
    virtual ~NodesOnBisectingPlanes();


    // This over-rides the base version so that the nodes are updated in the
    // order of dependence.
    virtual void SetNodeInAdjustmentOrder( size_t const adjustmentOrderIndex,
                                   std::vector< double > const& nodeAsVector )
    { pathNodes[ adjustmentOrder[ adjustmentOrderIndex ] ] = nodeAsVector; }


  protected:
    std::vector< size_t > adjustmentOrder;
    std::vector< std::pair< size_t, size_t > > sideNodeIndices;


    // This should add the perpendicular component from the parameterization
    // given by nodeParameterization along with startNode and endNode to
    // nodeVector.
    virtual void
    AddTransformedNode( std::vector< double > const& nodeVector,
                        std::vector< double > const& startNode,
                        std::vector< double > const& endNode,
                     std::vector< double > const& nodeParameterization ) const;

    // This should return the false-vacuum-side node of the pair of nodes
    // from which the node at adjustmentOrderIndex should be set.
    std::vector< double > const&
    FalseSideNode( size_t const adjustmentOrderIndex,
                   std::vector< std::vector< double > > const& nodeSet
                                                            = pathNodes ) const
    { return nodeSet[ sideNodeIndices[ adjustmentOrderIndex ].first ]; }

    // This should return the true-vacuum-side node of the pair of nodes
    // from which the node at nodeIndex should be set.
    std::vector< double > const&
    TrueSideNode( size_t const adjustmentOrderIndex,
                  std::vector< std::vector< double > > const& nodeSet
                                                       = pathNodes ) const
    { return nodeSet[ sideNodeIndices[ adjustmentOrderIndex].second ]; }

    // This should return the fraction along the node difference vector that
    // the rotated plane should be shifted appropriate for
    // pathNodes[ nodeIndex ].
    double ShiftFraction( size_t const nodeIndex ) const{ return 0.5; }
  };

} /* namespace VevaciousPlusPlus */
#endif /* NODESONBISECTINGPLANES_HPP_ */
