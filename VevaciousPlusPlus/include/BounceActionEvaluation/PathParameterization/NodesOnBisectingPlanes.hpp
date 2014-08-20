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
    // the appropriate element of reflectionMatrices.
    virtual void SetNodeFromNodeVector( size_t const nodeIndex,
                                   std::vector< double > const& nodeAsVector );


    // This returns the correct node index in path order corresponding to the
    // index in adjustment order so that the middle node is moved 1st, then the
    // nodes between the ends and the middle, and then the nodes in between
    // those, and so on.
    virtual size_t
    PathIndexFromAdjustmentIndex( size_t const adjustmentOrderIndex ) const
    { return adjustmentOrder[ adjustmentOrderIndex ]; }

    // This converts sideNodeIndices so that each node depends on its nearest
    // neighbors for its bisection, rather than the initial set-up. It also
    // updates reflectionMatrices appropriately.
    virtual void ConvertToTuningOrder();


  protected:
    std::vector< size_t > adjustmentOrder;
    std::vector< std::pair< size_t, size_t > > sideNodeIndices;
    std::vector< Eigen::MatrixXd > reflectionMatrices;
    // The nodes are parameterized as vectors perpendicular to the axis of the
    // 0th field, which are then transformed by the Householder reflection that
    // reflects the axis of the 0th field onto the vector of the difference
    // between the nodes


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

    // This sets reflectionMatrices[ nodeIndex ] to be a matrix that rotates the
    // vector difference from FalseSideNode( nodeIndex, pathNodes ) ) to
    // TrueSideNode( nodeIndex, pathNodes ) to lie along the axis of
    // referenceField.
    void UpdateRotationMatrix( size_t const nodeIndex );
  };



  // This sets pathNodes[ nodeIndex ] to be nodeAsVector, and also updates the
  // appropriate elements of reflectionMatrices.
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

  // This converts sideNodeIndices so that each node depends on its nearest
  // neighbors for its bisection, rather than the initial set-up. It also
  // updates reflectionMatrices appropriately.
  inline void NodesOnBisectingPlanes::ConvertToTuningOrder()
  {
    for( size_t sideIndex( 1 );
         sideIndex < ( pathNodes.size() - 1 );
         ++sideIndex )
    {
      sideNodeIndices[ sideIndex ].first = ( sideIndex - 1 );
      sideNodeIndices[ sideIndex ].second = ( sideIndex + 1 );
      UpdateRotationMatrix( sideIndex );
    }
  }

  // This ensures that the rotation matrices are set up.
  inline void NodesOnBisectingPlanes::FinishUpdatingForNewVacua()
  {
    UpdateRotationMatrix( PathIndexFromAdjustmentIndex( 0 ) );
    for( size_t matrixIndex( 0 );
         matrixIndex < reflectionMatrices.size();
         ++matrixIndex )
    {
      if( matrixIndex !=  PathIndexFromAdjustmentIndex( 0 ) )
      {
        reflectionMatrices[ matrixIndex ]
        = reflectionMatrices[  PathIndexFromAdjustmentIndex( 0 ) ];
      }
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "NodesOnBisectingPlanes::FinishUpdatingForNewVacua() finished."
    << " rotationMatrices = {"
    << std::endl;
    for( size_t matrixIndex( 0 );
         matrixIndex < reflectionMatrices.size();
         ++matrixIndex )
    {
      if( matrixIndex > 0 )
      {
        std::cout << "-------------" << std::endl;
      }
      std::cout << reflectionMatrices[ matrixIndex ] << std::endl;
    }
    std::cout << "}" << std::endl;
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESONBISECTINGPLANES_HPP_ */
