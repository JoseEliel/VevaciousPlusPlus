/*
 * NodesOnPlanes.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef NODESONPLANES_HPP_
#define NODESONPLANES_HPP_

#include "CommonIncludes.hpp"
#include "NodesFromParameterization.hpp"

namespace VevaciousPlusPlus
{

  class NodesOnPlanes : public NodesFromParameterization
  {
  public:
    NodesOnPlanes( std::vector< double > const& falseVacuum,
                   std::vector< double > const& trueVacuum,
                   size_t const numberOfIntermediateNodes );
    virtual
    ~NodesOnPlanes();


    // This returns the node at nodeIndex based on nodeParameterization, unless
    // it would mean moving the true vacuum or false vacuum, in which case an
    // exception is thrown. The node is set to be on a plane perpendicular
    // to two other nodes (decided by the derived class), by
    // nodeParameterization being taken to be a vector on a plane with the
    // field with index referenceField being zero, through
    // ProjectPerpendicularToAndShift.
    virtual std::vector< double > PathNode( size_t const nodeIndex,
                      std::vector< double > const& nodeParameterization ) const
    { return ProjectPerpendicularToAndShift( FalseSideNode( nodeIndex ),
                                             TrueSideNode( nodeIndex ),
                                             nodeParameterization,
                                             ShiftFraction( nodeIndex ) ); }


  protected:
    size_t referenceField;
    std::vector< double > const zeroParameterization;


    // This returns the vector sum of startNode plus shiftFraction times the
    // difference between startNode and endNode, plus nodeInPlane transformed
    // by TransformNodeInPlane to be perpendicular to the difference between
    // startNode and endNode.
    std::vector< double >
    ProjectPerpendicularToAndShift( std::vector< double > const& startNode,
                                    std::vector< double > const& endNode,
                                    std::vector< double > const& nodeInPlane,
                                    double const shiftFraction ) const;

    // This should return the false-vacuum-side node of the pair of nodes
    // from which the node at nodeIndex should be set.
    std::vector< double > const&
    FalseSideNode( size_t const nodeIndex ) const = 0;

    // This should return the true-vacuum-side node of the pair of nodes
    // from which the node at nodeIndex should be set.
    std::vector< double > const&
    TrueSideNode( size_t const nodeIndex ) const = 0;

    // This should return the fraction along the node difference vector that
    // the rotated plane should be shifted appropriate for
    // pathNodes[ nodeIndex ].
    size_t ShiftFraction( size_t const nodeIndex ) const = 0;
  };


  // This returns the node at nodeIndex based on nodeParameterization, unless
  // it would mean moving the true vacuum or false vacuum, in which case an
  // exception is thrown. The node is set to be on a plane perpendicular
  // to two other nodes (decided by the derived class), by
  // nodeParameterization being taken to be a vector on a plane with the
  // field with index referenceField being zero, through
  // ProjectPerpendicularToAndShift.
  virtual std::vector< double >
  NodesOnPlanes::PathNode( size_t const nodeIndex,
                      std::vector< double > const& nodeParameterization ) const
  {
    if( !( ( nodeIndex >= 1 )
           &&
           ( nodeIndex <= numberOfIntermediateNodes ) ) )
    {
      std::stringstream errorStream;
      errorStream << "Can only set nodes 1 to " << numberOfIntermediateNodes
      << ", not node " << nodeIndex << "!";
      throw std::out_of_range( errorStream.str() );
    }
    return ProjectPerpendicularToAndShift( FalseSideNode( nodeIndex ),
                                           TrueSideNode( nodeIndex ),
                                           nodeParameterization,
                                           ShiftFraction( nodeIndex ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESONPLANES_HPP_ */
