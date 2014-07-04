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


    // This return all the nodes based on the numbers given in
    // pathParameterization. The nodes returned should be ordered in the
    // sequence that they are visited in the path from the false vacuum to the
    // true vacuum, with front() being the false vacuum and back() being the
    // true vacuum.
    virtual std::vector< std::vector< double > >
    PathNodes( std::vector< double > const& pathParameterization ) const = 0;

    // This throws an exception by default, but can be over-ridden.
    virtual std::vector< double > PathNode( size_t const nodeIndex,
                      std::vector< double > const& nodeParameterization ) const
    { throw std::out_of_range( "No functionality to set just one node in"
                               " NodesFromParameterization base class!" ); }


  protected:
    size_t numberOfFields;
    size_t numberOfIntermediateNodes;
    std::vector< double > const& falseVacuum;
    std::vector< double > const& trueVacuum;


    // This should return the transformation of nodeInPlane as a vector with
    // numberOfFields components into a vector perpendicular to vector
    // difference (endNode - startNode).
    virtual std::vector< double >
    TransformNodeInPlane( std::vector< double > const& startNode,
                          std::vector< double > const& endNode,
                          std::vector< double > const& nodeInPlane ) const = 0;
  };




  // This sets pathNodes[ nodeIndex ] to be the vector sum of startNode plus
  // shiftFraction times the difference between startNode and endNode, plus
  // nodeInPlane rotated to be perpendicular to the difference between
  // startNode and endNode.
  inline std::vector< double > NodesOnPlanes::ProjectPerpendicularToAndShift(
                                        std::vector< double > const& startNode,
                                          std::vector< double > const& endNode,
                                      std::vector< double > const& nodeInPlane,
                                             double const shiftFraction ) const
  {
    std::vector< double > returnNode( TransformNodeInPlane( startNode,
                                                            endNode,
                                                            nodeInPlane ) );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      returnNode += ( ( 1.0 - shiftFraction ) * startNode[ fieldIndex ]
                      + ( shiftFraction * endNode[ fieldIndex ] ) );
    }
    return returnNode;
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESFROMPARAMETERIZATION_HPP_ */
