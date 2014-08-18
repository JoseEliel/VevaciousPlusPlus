/*
 * NodesOnParallelPlanes.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef NODESONPARALLELPLANES_HPP_
#define NODESONPARALLELPLANES_HPP_

#include "CommonIncludes.hpp"
#include "NodesOnPlanes.hpp"

namespace VevaciousPlusPlus
{

  class NodesOnParallelPlanes : public NodesOnPlanes
  {
  public:
    NodesOnParallelPlanes( size_t const numberOfFields,
                           size_t const numberOfIntermediateNodes );
    virtual ~NodesOnParallelPlanes();


  protected:
    // This takes nodeParameterization as a vector in the plane with field
    // referenceField = 0 and projects it onto the plane perpendicular to the
    // difference vector between the vacua, and adds that to nodeVector.
    virtual void AddTransformedNode( std::vector< double >& nodeVector,
                                     std::vector< double > const& startNode,
                                     std::vector< double > const& endNode,
                     std::vector< double > const& nodeParameterization ) const;

    // This returns the false vacuum node as the false-vacuum-side node of the
    // pair of nodes from which the node at adjustmentOrderIndex should be set.
    virtual std::vector< double > const&
    FalseSideNode( size_t const adjustmentOrderIndex,
                   std::vector< std::vector< double > > const& nodeSet ) const
    { return nodeSet.front(); }

    // This returns the true vacuum node as the true-vacuum-side node of the
    // pair of nodes from which the node at adjustmentOrderIndex should be set.
    virtual std::vector< double > const&
    TrueSideNode( size_t const adjustmentOrderIndex,
                  std::vector< std::vector< double > > const& nodeSet ) const
    { return nodeSet.back(); }

    // This returns the fraction along the difference vector between the vacua
    // that nodeIndex corresponds to.
    virtual double ShiftFraction( size_t const nodeIndex ) const
    { return ( static_cast<double>( nodeIndex )
               / static_cast<double>( numberOfIntermediateNodes + 1 ) ); }
  };

} /* namespace VevaciousPlusPlus */
#endif /* NODESONPARALLELPLANES_HPP_ */
