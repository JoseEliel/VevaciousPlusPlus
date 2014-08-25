/*
 * NodesFromParameterization.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef NODESFROMPARAMETERIZATION_HPP_
#define NODESFROMPARAMETERIZATION_HPP_

#include "CommonIncludes.hpp"
#include "Eigen/Dense"
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
    // nodeVector to be what the node at index nodeIndex would be for the
    // parameterization nodeParameterization. This throws an exception by
    // default.
    virtual void NodeProposal( std::vector< double >& nodeVector,
                               size_t const nodeIndex,
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

    // This returns the value of the index corresponding to the start of the
    // range of an index in adjustment order.
    virtual size_t AdjustmentOrderStartIndex() const{ return 0; }

    // This returns the value of the index corresponding to the end of the
    // range of an index in adjustment order.
    virtual size_t AdjustmentOrderEndIndex() const
    { return ( numberOfIntermediateNodes - 1 ); }

    // This returns false by default, but if a derived class has the
    // possibility of converting itself to a mode where the path can be further
    // refined, after say a coarse initial setup phase, then this should return
    // true.
    virtual bool HasFurtherRefinementMode() const{ return false; }

    // This does nothing by default, but is used by NodesOnBisectingPlanes to
    // change around the node dependences.
    virtual void ConvertToRefinementOrder(){}

    // This sets pathNodes[ nodeIndex ] to be nodeAsVector. (If the adjustment
    // order is different from the order the path visits the nodes in, the
    // calling code can account for that by using
    // PathIndexFromAdjustmentIndex( adjustmentOrderIndex ) first.)
    virtual void
    SetNodeFromNodeVector( size_t const nodeIndex,
                           std::vector< double > const& nodeAsVector )
    { pathNodes[ nodeIndex ] = nodeAsVector; }

    // This sets pathNodes[ nodeIndex ] to be the node parameterized by
    // nodeParameterization. (If the adjustment order is different from the
    // order the path visits the nodes in, the calling code can account for
    // that by using PathIndexFromAdjustmentIndex( adjustmentOrderIndex )
    // first.)
    void SetNodeFromParameterization( size_t const nodeIndex,
                           std::vector< double > const& nodeParameterization )
    { NodeProposal( pathNodes[ nodeIndex ],
                    nodeIndex,
                    nodeParameterization ); }

    std::vector< std::vector< double > >const&
    PathNodes() const{ return pathNodes; }

    size_t NumberOfFields() const{ return numberOfFields; }

    size_t NumberOfVaryingNodes() const{ return numberOfIntermediateNodes; }

    virtual size_t NumberOfPathNodes() const
    { return ( numberOfIntermediateNodes + 2 ); }

    std::vector< double > const& ZeroParameterization() const
    { return zeroFullParameterization; }

    std::vector< double > const& ZeroNodeParameterization() const
    { return zeroNodeParameterization; }

    std::vector< double > const& InitialStepSizes() const
    { return initialStepSizes; }

    // This just assumes that the nodes can be adjusted in the order which they
    // are visited in along the path (skipping the false vacuum at the front of
    // pathNodes), but can be over-ridden in derived classes.
    virtual size_t
    PathIndexFromAdjustmentIndex( size_t const adjustmentOrderIndex ) const
    { return ( adjustmentOrderIndex + 1 ); }

    // This just throws an exception if pathOrderIndex is out of the range of
    // pathNodes, or if it would be for the front or back of pathNodes as then
    // the vacua would be moved and that is not a good idea.
    void ValidPathIndex( size_t const pathOrderIndex ) const;

    // This returns true by default, but if a derived class is such that once a
    // node is set in adjustment order, it doesn't matter what nodes are set
    // afterwards, that set node won't need to be moved again based on where
    // the subsequent nodes ended up, then the derived class should over-ride
    // this to return false.
    virtual bool PreviousNodesDependOnSubsequentNodesInAdjustmentOrder() const
    { return true; }


  protected:
    size_t numberOfFields;
    size_t numberOfIntermediateNodes;
    std::vector< std::vector< double > > pathNodes;
    std::vector< double > zeroFullParameterization;
    std::vector< double > zeroNodeParameterization;
    std::vector< double > initialStepSizes;


    // This sets differenceVector[i] = ( endVector[i] - startVector[i] ) for
    // i running over all the indices of the vectors.
    void SetAsVectorDifference( std::vector< double >& differenceVector,
                                std::vector< double > const& startVector,
                                std::vector< double > const& endVector ) const;

    // This sets reflectionMatrix to be the Householder reflection which
    // reflects the axis of field referenceAxis to lie along targetVector.
    void SetAsHouseholderReflectionFromAxisToVector(
                                             Eigen::MatrixXd& reflectionMatrix,
                                                    size_t const referenceAxis,
                                   std::vector< double > const& targetVector );
  };




  // This resets the NodesFromParameterization so that it will produce
  // TunnelPath*s that parameterize the path between the given vacua.
  inline void
  NodesFromParameterization::SetVacua( PotentialMinimum const& falseVacuum,
                                       PotentialMinimum const& trueVacuum )
  {
    pathNodes.front() = falseVacuum.FieldConfiguration();
    pathNodes.back() = trueVacuum.FieldConfiguration();
    SetInitialParameterizationAndStepSizes( zeroFullParameterization,
                                            initialStepSizes );
  }

  // This just throws an exception if pathOrderIndex is out of the range of
  // pathNodes, or if it would be for the front or back of pathNodes as then
  // the vacua would be moved and that is not a good idea.
  inline void NodesFromParameterization::ValidPathIndex(
                                            size_t const pathOrderIndex ) const
  {
    if( !( ( pathOrderIndex > 0 )
           &&
           ( pathOrderIndex <= numberOfIntermediateNodes ) ) )
    {
      std::stringstream errorStream;
      errorStream << "Can only set (in \"path order\") nodes 1 to "
      << numberOfIntermediateNodes << ", not node " << pathOrderIndex << "!";
      throw std::out_of_range( errorStream.str() );
    }
  }


  // This sets differenceVector[i] = ( endVector[i] - startVector[i] ) for
  // i running over all the indices of the vectors.
  inline void NodesFromParameterization::SetAsVectorDifference(
                                       std::vector< double >& differenceVector,
                                      std::vector< double > const& startVector,
                                 std::vector< double > const& endVector ) const
  {
    differenceVector = endVector;
    for( size_t fieldIndex( 0 );
         fieldIndex < startVector.size();
         ++fieldIndex )
    {
      differenceVector[ fieldIndex ] -= startVector[ fieldIndex ];
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESFROMPARAMETERIZATION_HPP_ */
