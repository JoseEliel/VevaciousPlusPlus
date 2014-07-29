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
    virtual ~NodesOnPlanes();

    // This puts all the nodes based on the numbers given in
    // pathParameterization into pathNodes, ordered in the sequence that they
    // are visited in the path from the false vacuum to the true vacuum, with
    // pathNodes.front() being the false vacuum and pathNodes.back() being the
    // true vacuum. (It would be nice to just return a new vector, but we're
    // avoiding relying on C++11-compliant compilers.)
    void PathNodeSet( std::vector< std::vector< double > >& pathNodes,
                     std::vector< double > const& pathParameterization ) const;

    // This sets pathParameterization to be the set of zero nodes in planes,
    // and initialStepSizes to be the same size, with each element being half
    // the Euclidean length of the difference between the vacua in field space.
    virtual void SetInitialParameterizationAndStepSizes(
                                   std::vector< double >& pathParameterization,
                               std::vector< double >& initialStepSizes ) const;

    // This sets nodeVector to be the vector which would represent the node at
    // pathNodes[ adjustmentOrderIndex ] if it were to be set based on
    // nodeParameterization and the rest of the nodes in pathNodes, unless it
    // would mean moving the true vacuum or false vacuum, in which case an
    // exception is thrown. The node is set to be on a plane perpendicular to
    // two other nodes (decided by the derived class), by
    // nodeParameterization being taken to be a vector on a plane with the
    // field with index referenceField being zero, through
    // ProjectPerpendicularToAndShift.
    virtual void NodeProposal( std::vector< double >& nodeVector,
                               size_t const adjustmentOrderIndex,
                      std::vector< double > const& nodeParameterization ) const
    { TransformPerpendicularToAndShift( nodeVector,
                                        FalseSideNode( adjustmentOrderIndex,
                                                       pathNodes ),
                                        TrueSideNode( adjustmentOrderIndex,
                                                      pathNodes ),
                                        nodeParameterization,
                                     ShiftFraction( adjustmentOrderIndex ) ); }

    // This sets nodeParameterization to be the sub-vector of
    // pathParameterization which parameterizes the node at index nodeIndex.
    virtual void ExtractSingleNodeParameterization(
                                   std::vector< double >& nodeParameterization,
                                                    size_t const nodeIndex,
                     std::vector< double > const& pathParameterization ) const
    { nodeParameterization.assign( ( pathParameterization.begin()
                                 + ( nodeIndex * numberOfParametersPerNode ) ),
                                   ( pathParameterization.begin()
                     + ( ( nodeIndex + 1 ) * numberOfParametersPerNode ) ) ); }

    // This sets nodeParameterization to be appropriate for a set of initial
    // step sizes for Minuit2. This throws an exception by default.
    virtual void SetSingleNodeStepSizes(
                               std::vector< double >& initialStepSizes ) const;


  protected:
    size_t referenceField;
    size_t numberOfParametersPerNode;


    // This should add the perpendicular component from the parameterization
    // given by nodeParameterization along with startNode and endNode to
    // nodeVector.
    virtual void
    AddTransformedNode( std::vector< double >& nodeVector,
                        std::vector< double > const& startNode,
                        std::vector< double > const& endNode,
                 std::vector< double > const& nodeParameterization ) const = 0;

    // This sets nodeVector to be the vector sum of startNode plus
    // shiftFraction times the difference between startNode and endNode, plus
    // the node given by nodeParameterization transformed to be perpendicular
    // to the difference between startNode and endNode.
    void TransformPerpendicularToAndShift( std::vector< double >& nodeVector,
                                        std::vector< double > const& startNode,
                                          std::vector< double > const& endNode,
                             std::vector< double > const& nodeParameterization,
                                           double const shiftFraction ) const;

    // This should return the false-vacuum-side node of the pair of nodes
    // from which the node at adjustmentOrderIndex should be set.
    virtual std::vector< double > const&
    FalseSideNode( size_t const adjustmentOrderIndex,
               std::vector< std::vector< double > > const& nodeSet ) const = 0;

    // This should return the true-vacuum-side node of the pair of nodes
    // from which the node at nodeIndex should be set.
    virtual std::vector< double > const&
    TrueSideNode( size_t const adjustmentOrderIndex,
               std::vector< std::vector< double > > const& nodeSet ) const = 0;

    // This should return the fraction along the node difference vector that
    // the rotated plane should be shifted appropriate for
    // pathNodes[ nodeIndex ].
    virtual double ShiftFraction( size_t const nodeIndex ) const = 0;
  };




  // This puts all the nodes based on the numbers given in
  // pathParameterization into pathNodes, ordered in the sequence that they
  // are visited in the path from the false vacuum to the true vacuum, with
  // pathNodes.front() being the false vacuum and pathNodes.back() being the
  // true vacuum. (It would be nice to just return a new vector, but we're
  // avoiding relying on C++11-compliant compilers.)
  void
  NodesOnPlanes::PathNodeSet( std::vector< std::vector< double > >& pathNodes,
                      std::vector< double > const& pathParameterization ) const
  {
    pathNodes = this->pathNodes;
    std::vector< double > nodeParameterization( numberOfParametersPerNode );
    std::vector< double >::const_iterator
    currentParameterizationStart( pathParameterization.begin() );
    for( size_t nodeIndex( 1 );
         nodeIndex <= numberOfIntermediateNodes;
         ++nodeIndex )
    {
      nodeParameterization.assign( currentParameterizationStart,
                ( currentParameterizationStart + numberOfParametersPerNode ) );
      TransformPerpendicularToAndShift( pathNodes[ nodeIndex ],
                                        FalseSideNode( nodeIndex,
                                                       pathNodes ),
                                        TrueSideNode( nodeIndex,
                                                      pathNodes ),
                                        nodeParameterization,
                                        ShiftFraction( nodeIndex ) );
      currentParameterizationStart += numberOfParametersPerNode;
    }
  }

  // This sets pathParameterization to be the set of zero nodes in planes,
  // and initialStepSizes to be the same size, with each element being half
  // the Euclidean length of the difference between the vacua in field space.
  inline void NodesOnPlanes::SetInitialParameterizationAndStepSizes(
                                  std::vector< double >& pathParameterization,
                                std::vector< double >& initialStepSizes ) const
  {
    pathParameterization.assign( ( numberOfIntermediateNodes
                                   * numberOfParametersPerNode ),
                                 0.0 );
    double fieldLengthSquared( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double const fieldDifference( pathNodes.back()[ fieldIndex ]
                                    - pathNodes.front()[ fieldIndex ] );
      fieldLengthSquared += ( fieldDifference * fieldDifference );
    }
    initialStepSizes.assign( ( numberOfIntermediateNodes
                               * numberOfParametersPerNode ),
                              ( 0.5 * sqrt( fieldLengthSquared ) ) );
  }

  // This sets nodeParameterization to be numberOfParametersPerNode elements,
  // each being half the Euclidean length of the difference between the vacua
  // in field space.
  inline void NodesOnPlanes::SetSingleNodeStepSizes(
                                std::vector< double >& initialStepSizes ) const
  {
    double fieldLengthSquared( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double const fieldDifference( pathNodes.back()[ fieldIndex ]
                                    - pathNodes.front()[ fieldIndex ] );
      fieldLengthSquared += ( fieldDifference * fieldDifference );
    }
    initialStepSizes.assign( numberOfParametersPerNode,
                             ( 0.5 * sqrt( fieldLengthSquared ) ) );
  }

  // This sets nodeVector to be the vector sum of startNode plus
  // shiftFraction times the difference between startNode and endNode, plus
  // the node given by nodeParameterization transformed to be perpendicular
  // to the difference between startNode and endNode.
  inline void NodesOnPlanes::TransformPerpendicularToAndShift(
                                             std::vector< double >& nodeVector,
                                        std::vector< double > const& startNode,
                                          std::vector< double > const& endNode,
                             std::vector< double > const& nodeParameterization,
                                             double const shiftFraction ) const
  {
    // We actually take the fraction along the straight line between nodes (the
    // "shift") before adding the perpendicular component.
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      nodeVector[ fieldIndex ]
      = ( ( 1.0 - shiftFraction ) * startNode[ fieldIndex ]
          + ( shiftFraction * endNode[ fieldIndex ] ) );
    }
    if( nodeParameterization != zeroParameterization )
    {
      AddTransformedNode( nodeVector,
                          startNode,
                          endNode,
                          nodeParameterization );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESONPLANES_HPP_ */
