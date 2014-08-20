/*
 * NodesOnParallelPlanes.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef NODESONPARALLELPLANES_HPP_
#define NODESONPARALLELPLANES_HPP_

#include "CommonIncludes.hpp"
#include "Eigen/Dense"
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
    Eigen::MatrixXd reflectionMatrix;


    // This takes nodeParameterization as a vector in the plane with the 0th
    // field being 0.0 and projects it onto the plane perpendicular to the
    // difference vector between the vacua, and adds that to nodeVector.
    virtual void AddTransformedNode( std::vector< double >& nodeVector,
                                     size_t const nodeIndex,
                     std::vector< double > const& nodeParameterization ) const;

    // This returns the false vacuum node as the false-vacuum-side node of the
    // pair of nodes from which the node at nodeIndex should be set.
    virtual std::vector< double > const&
    FalseSideNode( size_t const nodeIndex,
                   std::vector< std::vector< double > > const& nodeSet ) const
    { return nodeSet.front(); }

    // This returns the true vacuum node as the true-vacuum-side node of the
    // pair of nodes from which the node at nodeIndex should be set.
    virtual std::vector< double > const&
    TrueSideNode( size_t const nodeIndex,
                  std::vector< std::vector< double > > const& nodeSet ) const
    { return nodeSet.back(); }

    // This returns the fraction along the difference vector between the vacua
    // that nodeIndex corresponds to.
    virtual double ShiftFraction( size_t const nodeIndex ) const
    { return ( static_cast< double >( nodeIndex )
               / static_cast< double >( numberOfIntermediateNodes + 1 ) ); }

    // This ensures that the reflection matrix is set up.
    virtual void FinishUpdatingForNewVacua();


    // This is the old version of AddTransformedNode, before moving to using
    // the Householder reflection matrix.
    // This takes nodeParameterization as a vector in the plane with field
    // referenceField = 0 and projects it onto the plane perpendicular to the
    // difference vector between the vacua, and adds that to nodeVector.
    virtual void AddProjectedNode( std::vector< double >& nodeVector,
                                   size_t const nodeIndex,
                     std::vector< double > const& nodeParameterization ) const;
  };




  // This ensures that the reflection matrix is set up.
  inline void NodesOnParallelPlanes::FinishUpdatingForNewVacua()
  {
    std::vector< double > vacuumDifference;
    SetAsVectorDifference( vacuumDifference,
                           pathNodes.front(),
                           pathNodes.back() );
    SetAsHouseholderReflectionFromAxisToVector( reflectionMatrix,
                                                0,
                                                vacuumDifference );
    // For convenience, the 0th field is taken to be the 0.0 in the
    // parameterization. Taking anything else would complicate the rest of the
    // functions.
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESONPARALLELPLANES_HPP_ */
