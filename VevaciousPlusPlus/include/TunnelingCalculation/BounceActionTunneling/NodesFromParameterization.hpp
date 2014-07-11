/*
 * NodesFromParameterization.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef NODESFROMPARAMETERIZATION_HPP_
#define NODESFROMPARAMETERIZATION_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class NodesFromParameterization
  {
  public:
    NodesFromParameterization( std::vector< double > const& falseVacuum,
                               std::vector< double > const& trueVacuum,
                               size_t const numberOfIntermediateNodes );
    virtual ~NodesFromParameterization();


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

    // This sets pathNodes[ adjustmentOrderIndex ] to be nodeAsVector as long
    // as adjustmentOrderIndex is >= 1 and <= numberOfIntermediateNodes so that
    // it doesn't try to overwrite the vacua or anything out of the range of
    // pathNodes, throwing an out-of-range exception otherwise. It can be
    // over-ridden if necessary.
    virtual void SetNodeInAdjustmentOrder( size_t const adjustmentOrderIndex,
                                   std::vector< double > const& nodeAsVector );

    size_t NumberOfFields() const{ return numberOfFields; }

    size_t NumberOfVaryingNodes() const{ return numberOfIntermediateNodes; }

    virtual size_t NumberOfPathNodes() const
    { return ( numberOfIntermediateNodes + 2 ); }

    std::vector< double > const& ZeroParameterization() const
    { return zeroParameterization; }


  protected:
    size_t numberOfFields;
    size_t numberOfIntermediateNodes;
    std::vector< std::vector< double > > pathNodes;
    std::vector< double > zeroParameterization;
  };




  // This sets pathNodes[ adjustmentOrderIndex ] to be nodeAsVector as long
  // as adjustmentOrderIndex is >= 1 and <= numberOfIntermediateNodes so that
  // it doesn't try to overwrite the vacua or anything out of the range of
  // pathNodes, throwing an out-of-range exception otherwise. It can be
  // over-ridden if necessary.
  inline void
  NodesOnPlanes::SetNodeInAdjustmentOrder( size_t const adjustmentOrderIndex,
                                    std::vector< double > const& nodeAsVector )
  {
    if( ( adjustmentOrderIndex > 0 )
        &&
        ( adjustmentOrderIndex < ( pathNodes.size() - 1 ) ) )
    {
      pathNodes[ adjustmentOrderIndex ] = nodeAsVector;
    }
    else
    {
      std::stringstream errorStream;
      errorStream << "Can only set nodes 1 to " << numberOfIntermediateNodes
      << ", not node " << adjustmentOrderIndex << "!";
      throw std::out_of_range( errorStream.str() );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESFROMPARAMETERIZATION_HPP_ */
