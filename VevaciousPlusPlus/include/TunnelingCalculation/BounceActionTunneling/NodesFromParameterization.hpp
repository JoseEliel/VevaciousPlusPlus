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


    // This should return all the nodes based on the numbers given in
    // pathParameterization. The nodes returned should be ordered in the
    // sequence that they are visited in the path from the false vacuum to the
    // true vacuum, with front() being the false vacuum and back() being the
    // true vacuum.
    virtual std::vector< std::vector< double > >
    PathNodeSet( std::vector< double > const& pathParameterization ) const = 0;

    // This throws an exception by default, but can be over-ridden.
    virtual std::vector< double >
    NodeProposal( size_t const adjustmentOrderIndex,
                  std::vector< double > const& nodeParameterization ) const
    { throw std::out_of_range( "No functionality to set just one node in"
                               " NodesFromParameterization base class!" ); }

    // This sets pathNodes[ adjustmentOrderIndex ] to be nodeAsVector as long
    // as adjustmentOrderIndex is >= 1 and <= numberOfIntermediateNodes so that
    // it doesn't try to overwrite the vacua or anything out of the range of
    // pathNodes, throwing an out-of-range exception otherwise. It can be
    // over-ridden if necessary.
    virtual void SetNodeInAdjustmentOrder( size_t const adjustmentOrderIndex,
                                   std::vector< double > const& nodeAsVector );


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
