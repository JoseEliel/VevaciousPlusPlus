/*
 * NodesFromParameterization.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef NODESFROMPARAMETERIZATION_HPP_
#define NODESFROMPARAMETERIZATION_HPP_

#include "CommonIncludes.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"

namespace VevaciousPlusPlus
{

  class NodesFromParameterization
  {
  public:
    NodesFromParameterization( size_t const numberOfFields,
                               size_t const numberOfIntermediateNodes );
    virtual ~NodesFromParameterization();


    // This resets the NodesFromParameterization so that it will produce
    // TunnelPath*s that parameterize the path between the given vacua.
    virtual void SetVacua( PotentialMinimum const& falseVacuum,
                           PotentialMinimum const& trueVacuum );

    // This should put all the nodes based on the numbers given in
    // pathParameterization into pathNodes, ordered in the sequence that they
    // are visited in the path from the false vacuum to the true vacuum, with
    // pathNodes.front() being the false vacuum and pathNodes.back() being the
    // true vacuum. (It would be nice to just return a new vector, but we're
    // avoiding relying on C++11-compliant compilers.)
    virtual void PathNodeSet( std::vector< std::vector< double > >& pathNodes,
                 std::vector< double > const& pathParameterization ) const = 0;

    // This should set pathParameterization to be the default initial
    // parameterization for the derived class, and initialStepSizes to be the
    // default initial step sizes that would suit Minuit2.
    virtual void SetInitialParameterizationAndStepSizes(
                                   std::vector< double >& pathParameterization,
                           std::vector< double >& initialStepSizes ) const = 0;

    // This could construct a parameterization from the given nodes if
    // necessary, but right now it's not necessary, so passing in an empty
    // vector is sufficient for the cases where the TunnelPath is being created
    // directly from a set of nodes rather than a parameterization setting up a
    // set of nodes which then creates the TunnelPath.
    virtual std::vector< double > ParameterizationForNodes(
                  std::vector< std::vector< double > > const& pathNodes ) const
    { return std::vector< double >( 0 ); }

    // If the derived class can vary one node at a time, this should set
    // nodeVector to be what the node at index adjustmentOrderIndex would be
    // for the parameterization nodeParameterization. This throws an exception
    // by default.
    virtual void NodeProposal( std::vector< double >& nodeVector,
                               size_t const adjustmentOrderIndex,
                      std::vector< double > const& nodeParameterization ) const
    { throw std::out_of_range( "No functionality to set just one node in"
                               " NodesFromParameterization base class!" ); }

    // If the derived class can vary one node at a time, this should set
    // nodeParameterization to be the sub-vector of pathParameterization which
    // parameterizes the node at index nodeIndex. This throws an exception by
    // default.
    virtual void ExtractSingleNodeParameterization(
                                   std::vector< double >& nodeParameterization,
                                                    size_t const nodeIndex,
                      std::vector< double > const& pathParameterization ) const
    { throw std::out_of_range( "No functionality to set just one node in"
                               " NodesFromParameterization base class!" ); }

    // If the derived class can vary one node at a time, this should set
    // nodeParameterization to be appropriate for a set of initial step sizes
    // for Minuit2. This throws an exception by default.
    virtual void
    SetSingleNodeStepSizes( std::vector< double >& initialStepSizes ) const
    { throw std::out_of_range( "No functionality to set just one node in"
                               " NodesFromParameterization base class!" ); }

    // This sets
    // pathNodes[ PathIndexFromAdjustmentIndex( adjustmentOrderIndex ) ] ] to
    // be nodeAsVector as long as the index is valid. (If the adjustment order
    // is different from the order the path visits the nodes in,
    // PathIndexFromAdjustmentIndex can account for that.)
    virtual void
    SetNodeInAdjustmentOrderFromNodeVector( size_t const adjustmentOrderIndex,
                                    std::vector< double > const& nodeAsVector )
    { pathNodes[ PathIndexFromAdjustmentIndex( adjustmentOrderIndex ) ]
      = nodeAsVector; }

    // This sets
    // pathNodes[ PathIndexFromAdjustmentIndex( adjustmentOrderIndex ) ] ] to
    // be the node parameterized by nodeParameterization as long as the index
    // is valid. (If the adjustment order is different from the order the path
    // visits the nodes in, PathIndexFromAdjustmentIndex can account for that.)
    void SetNodeInAdjustmentOrderFromParameterization(
                                             size_t const adjustmentOrderIndex,
                           std::vector< double > const& nodeParameterization )
    { NodeProposal( pathNodes[ PathIndexFromAdjustmentIndex(
                                                      adjustmentOrderIndex ) ],
                    adjustmentOrderIndex,
                    nodeParameterization ); }

    std::vector< std::vector< double > >const&
    PathNodes() const{ return pathNodes; }

    size_t NumberOfFields() const{ return numberOfFields; }

    size_t NumberOfVaryingNodes() const{ return numberOfIntermediateNodes; }

    virtual size_t NumberOfPathNodes() const
    { return ( numberOfIntermediateNodes + 2 ); }

    std::vector< double > const& ZeroParameterization() const
    { return zeroParameterization; }

    std::vector< double > const& InitialStepSizes() const
    { return initialStepSizes; }


  protected:
    size_t numberOfFields;
    size_t numberOfIntermediateNodes;
    std::vector< std::vector< double > > pathNodes;
    std::vector< double > zeroParameterization;
    std::vector< double > initialStepSizes;


    // This just assumes that the nodes can be adjusted in the order which they
    // are visited in along the path, but can be over-ridden in derived
    // classes. (It also calls ValidAdjustmentIndex before returning the
    // index.)
    virtual size_t
    PathIndexFromAdjustmentIndex( size_t const adjustmentOrderIndex ) const;

    // This just throws an exception if adjustmentOrderIndex is out of the
    // range of pathNodes, or if it would be for the front or back of pathNodes
    // as then the vacua would be moved and that is not a good idea.
    void ValidAdjustmentIndex( size_t const adjustmentOrderIndex ) const;
  };




  // This resets the NodesFromParameterization so that it will produce
  // TunnelPath*s that parameterize the path between the given vacua.
  inline void
  NodesFromParameterization::SetVacua( PotentialMinimum const& falseVacuum,
                                       PotentialMinimum const& trueVacuum )
  {
    pathNodes.front() = falseVacuum.FieldConfiguration();
    pathNodes.back() = trueVacuum.FieldConfiguration();
    SetInitialParameterizationAndStepSizes( zeroParameterization,
                                            initialStepSizes );
  }

  // This just assumes that the nodes can be adjusted in the order which they
  // are visited in along the path, but can be over-ridden in derived
  // classes.
  inline size_t NodesFromParameterization::PathIndexFromAdjustmentIndex(
                                      size_t const adjustmentOrderIndex ) const
  {
    ValidAdjustmentIndex( adjustmentOrderIndex );
    return adjustmentOrderIndex;
  }

  // This just throws an exception if adjustmentOrderIndex is out of the
  // range of pathNodes, or if it would be for the front or back of pathNodes
  // as then the vacua would be moved and that is not a good idea.
  inline void NodesFromParameterization::ValidAdjustmentIndex(
                                      size_t const adjustmentOrderIndex ) const
  {
    if( !( ( adjustmentOrderIndex > 0 )
           &&
           ( adjustmentOrderIndex < ( pathNodes.size() - 1 ) ) ) )
    {
      std::stringstream errorStream;
      errorStream << "Can only set nodes 1 to " << numberOfIntermediateNodes
      << ", not node " << adjustmentOrderIndex << "!";
      throw std::out_of_range( errorStream.str() );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESFROMPARAMETERIZATION_HPP_ */
