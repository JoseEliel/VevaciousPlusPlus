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
    NodesOnBisectingPlanes( size_t const numberOfFields,
                            size_t const numberOfIntermediateNodes );
    virtual ~NodesOnBisectingPlanes();


    // This sets pathNodes[ nodeIndex ] to be nodeAsVector, and also updates
    // the appropriate element of rotationMatrices.
    virtual void SetNodeFromNodeVector( size_t const nodeIndex,
                                   std::vector< double > const& nodeAsVector );


    // This returns the correct node index in path order corresponding to the
    // index in adjustment order so that the middle node is moved 1st, then the
    // nodes between the ends and the middle, and then the nodes in between
    // those, and so on.
    virtual size_t
    PathIndexFromAdjustmentIndex( size_t const adjustmentOrderIndex ) const
    { return adjustmentOrder[ adjustmentOrderIndex ]; }


  protected:
    std::vector< size_t > adjustmentOrder;
    std::vector< std::pair< size_t, size_t > > sideNodeIndices;
    std::vector< Eigen::MatrixXd > rotationMatrices;


    // This adds the perpendicular component from the parameterization given by
    // nodeParameterization along with nodeIndex to nodeVector.
    virtual void AddTransformedNode( std::vector< double >& nodeVector,
                                     size_t const nodeIndex,
                     std::vector< double > const& nodeParameterization ) const;

    // This returns the false-vacuum-side node of the pair of nodes from which
    // the node at nodeIndex should be set.
    virtual std::vector< double > const&
    FalseSideNode( size_t const nodeIndex,
                   std::vector< std::vector< double > > const& nodeSet ) const
    { return nodeSet[ sideNodeIndices[ nodeIndex ].first ]; }

    // This returns the true-vacuum-side node of the pair of nodes from which
    // the node at nodeIndex should be set.
    virtual std::vector< double > const&
    TrueSideNode( size_t const nodeIndex,
                  std::vector< std::vector< double > > const& nodeSet ) const
    { return nodeSet[ sideNodeIndices[ nodeIndex].second ]; }

    // This should return the fraction along the node difference vector that
    // the rotated plane should be shifted appropriate for
    // pathNodes[ nodeIndex ].
    virtual double ShiftFraction( size_t const nodeIndex ) const{ return 0.5; }

    // This ensures that the rotation matrices are set up.
    virtual void FinishUpdatingForNewVacua();

    // This sets rotationMatrices[ nodeIndex ] to be a matrix that rotates the
    // vector difference from FalseSideNode( nodeIndex, pathNodes ) ) to
    // TrueSideNode( nodeIndex, pathNodes ) to lie along the axis of
    // referenceField.
    void UpdateRotationMatrix( size_t const nodeIndex );
  };



  // This sets pathNodes[ nodeIndex ] to be nodeAsVector, and also updates the
  // appropriate elements of rotationMatrices.
  inline void
  NodesOnBisectingPlanes::SetNodeFromNodeVector( size_t const nodeIndex,
                                    std::vector< double > const& nodeAsVector )
  {
    pathNodes[ nodeIndex ] = nodeAsVector;
    for( size_t sideIndex( 1 );
         sideIndex < ( pathNodes.size() - 1 );
         ++sideIndex )
    {
      if( ( sideNodeIndices[ sideIndex ].first == nodeIndex )
          ||
          ( sideNodeIndices[ sideIndex ].second == nodeIndex ) )
      {
        UpdateRotationMatrix( sideIndex );
      }
    }
  }

  // This ensures that the rotation matrices are set up.
  inline void NodesOnBisectingPlanes::FinishUpdatingForNewVacua()
  {
    UpdateRotationMatrix( PathIndexFromAdjustmentIndex( 0 ) );
    for( size_t matrixIndex( 0 );
         matrixIndex < rotationMatrices.size();
         ++matrixIndex )
    {
      if( matrixIndex !=  PathIndexFromAdjustmentIndex( 0 ) )
      {
        rotationMatrices[ matrixIndex ]
        = rotationMatrices[  PathIndexFromAdjustmentIndex( 0 ) ];
      }
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESONBISECTINGPLANES_HPP_ */
