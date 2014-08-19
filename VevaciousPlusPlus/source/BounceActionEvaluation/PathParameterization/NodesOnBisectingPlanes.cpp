/*
 * NodesOnBisectingPlanes.cpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/NodesOnBisectingPlanes.hpp"

namespace VevaciousPlusPlus
{

  NodesOnBisectingPlanes::NodesOnBisectingPlanes( size_t const numberOfFields,
                                     size_t const numberOfIntermediateNodes ) :
    NodesOnPlanes( numberOfFields,
                   numberOfIntermediateNodes ),
    adjustmentOrder(),
    sideNodeIndices(),
    rotationMatrices()
  {
    // We ensure that the number of intermediate nodes is one less than an
    // integer power of two, using the argument numberOfIntermediateNodes as a
    // lower bound.
    this->numberOfIntermediateNodes = 1;
    size_t numberOfSplits( 0 );
    while( this->numberOfIntermediateNodes <= numberOfIntermediateNodes )
    {
      this->numberOfIntermediateNodes *= 2;
      ++numberOfSplits;
    }
    --(this->numberOfIntermediateNodes);

    // We ensure that pathNodes and rotationMatrices have the correct size.
    pathNodes.resize( this->numberOfIntermediateNodes + 2 );
    sideNodeIndices.resize( pathNodes.size() );
    rotationMatrices.resize( this->numberOfIntermediateNodes );
    // There are only rotation matrices for the variable nodes.

    // We need to set up the order in which the nodes will be set. The middle
    // node is set first, based on the vacua at the ends of pathNodes, then
    // the nodes half-way between each pair from the previous round are set,
    // and so on.
    // For example, seven varying nodes, with pathNodes[ 0 ] being the false
    // vacuum and pathNodes[ 8 ] being the true vacuum:
    // pathNodes[ 4 = 8/2 ] is set first, based on pathNodes[ 0 = 4 - 4 ]
    // and pathNodes[ 8 = 4 + 4 ], then pathNodes[ 2 = 4/2 * ( 1 + 2 * 0 ) ]
    // based on pathNodes[ 0 = 2 - 2 ] and pathNodes[ 4 = 2 + 2 ], and
    // pathNodes[ 6 = 4/2 * ( 1 + 2 * 1 ) ] based on pathNodes[ 4 = 6 - 2 ] and
    // pathNodes[ 8 = 6 + 2 ], then the last round is
    // pathNodes[ 1 = 2/2 * ( 1 + 2 * 0 ) ] based on pathNodes[ 0 = 1 - 1 ] and
    // pathNodes[ 2 = 1 + 1 ], pathNodes[ 3 = 2/2 * ( 1 + 2 * 1 ) ] based on
    // pathNodes[ 2 = 3 - 1 ] and pathNodes[ 4 = 3 + 1 ],
    // pathNodes[ 5 = 2/2 * ( 1 + 2 * 2 ) ] based on pathNodes[ 4 = 5 - 1 ] and
    // pathNodes[ 6 = 5 + 1 ], pathNodes[ 7 = 2/2 * ( 1 + 2 * 3 ) ] based on
    // pathNodes[ 6 = 7 - 1 ] and pathNodes[ 8 = 7 + 1 ].
    size_t segmentSize( pathNodes.size() - 1 );
    size_t newSegmentsInSplit( 1 );
    for( size_t splitCount( 0 );
         splitCount < numberOfSplits;
         ++splitCount )
    {
      segmentSize /= 2;
      for( size_t whichSegment( 0 );
           whichSegment < newSegmentsInSplit;
           ++whichSegment )
      {
        size_t const
        currentIndex( segmentSize * ( 1 + ( 2 * whichSegment ) ) );
        adjustmentOrder.push_back( currentIndex );
        sideNodeIndices[ currentIndex ].first = ( currentIndex - segmentSize );
        sideNodeIndices[ currentIndex ].second
        = ( currentIndex + segmentSize );
      }
      newSegmentsInSplit *= 2;
    }
  }

  NodesOnBisectingPlanes::~NodesOnBisectingPlanes()
  {
    // This does nothing.
  }


  // This adds the perpendicular component from the parameterization given by
  // nodeParameterization along with nodeIndex to nodeVector.
  void NodesOnBisectingPlanes::AddTransformedNode(
                                             std::vector< double >& nodeVector,
                                                   size_t const nodeIndex,
                      std::vector< double > const& nodeParameterization ) const
  {
    // The process is to create a vector with numberOfFields components out of
    // nodeParameterization which is in the appropriate plane, then rotate it
    // by the rotation which we use consistently to rotate the referenceField
    // axis to align with (endNode - startNode). There is not a unique rotation
    // if numberOfFields is larger than two, so we choose the easiest thing to
    // implement and keep it consistent.

    Eigen::VectorXd nodeInPlane( numberOfFields );
    nodeInPlane( referenceField ) = 0.0;
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex < referenceField )
      {
        nodeInPlane( fieldIndex ) = nodeParameterization[ fieldIndex ];
      }
      else if( fieldIndex > referenceField )
      {
        nodeInPlane( fieldIndex ) = nodeParameterization[ fieldIndex - 1 ];
      }
    }

    Eigen::VectorXd rotatedNode( rotationMatrices[ nodeIndex ] * nodeInPlane );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      nodeVector[ fieldIndex ] += rotatedNode( fieldIndex );
    }
  }

  // This sets rotationMatrices[ nodeIndex ] to be a matrix that rotates the
  // vector difference from FalseSideNode( nodeIndex, pathNodes ) ) to
  // TrueSideNode( nodeIndex, pathNodes ) to lie along the axis of
  // referenceField.
  void NodesOnBisectingPlanes::UpdateRotationMatrix( size_t const nodeIndex )
  {
    std::vector< double > const& startNode( FalseSideNode( nodeIndex,
                                                           pathNodes ) );
    std::vector< double > const& endNode( TrueSideNode( nodeIndex,
                                                        pathNodes ) );
    Eigen::MatrixXd& rotationMatrix( rotationMatrices[ nodeIndex ] );
    double firstDifference( endNode.front() - startNode.back() );
    for( size_t columnIndex( 0 );
         columnIndex < numberOfFields;
         ++columnIndex )
    {
      rotationMatrix( 0,
                      columnIndex ) = firstDifference;
    }
    for( size_t rowIndex( 1 );
         rowIndex < numberOfFields;
         ++rowIndex )
    {
      rotationMatrix( rowIndex,
                      referenceField )
      = ( endNode[ rowIndex ] - startNode[ rowIndex ] );
    }
    double const firstElementSquared( firstDifference * firstDifference );
    // We keep numeratorSquared as a negative number.
    double numeratorSquared( 0.0 );
    for( size_t columnIndex( 0 );
         columnIndex < ( numberOfFields - 1 );
         ++columnIndex )
    {
      size_t columnInsertionIndex( columnIndex );
      if( columnIndex >= referenceField )
      {
        ++columnInsertionIndex;
      }
      numeratorSquared -= ( rotationMatrix( columnIndex,
                                            referenceField )
                            * rotationMatrix( columnIndex,
                                              referenceField ) );
      for( size_t rowIndex( 1 );
           rowIndex <= columnIndex;
           ++rowIndex )
      {
        rotationMatrix( rowIndex,
                        columnInsertionIndex ) = rotationMatrix( rowIndex,
                                                              referenceField );
      }
      rotationMatrix( ( columnIndex + 1 ),
                      columnInsertionIndex )
      = ( numeratorSquared / rotationMatrix( ( columnIndex + 1 ),
                                             referenceField ) );
      for( size_t rowIndex( columnIndex + 2 );
           rowIndex < numberOfFields;
           ++rowIndex )
      {
        rotationMatrix( rowIndex,
                        columnInsertionIndex ) = 0.0;
      }
    }
    // At this point, the rows make a set of orthogonal vectors, but they are
    // not normalized to unit length. Currently -numeratorSquared is the length
    // of the first row, which is convenient.
    double inverseLengthSquared( -1.0 / numeratorSquared );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfFields;
         ++rowIndex )
    {
      rotationMatrix( rowIndex,
                      referenceField ) *= inverseLengthSquared;
    }
    double squaredSum( 0.0 );
    double nextSquare( firstElementSquared );
    for( size_t columnIndex( 0 );
         columnIndex < ( numberOfFields - 1 );
         ++columnIndex )
    {
      size_t columnNormalizationIndex( columnIndex );
      if( columnIndex >= referenceField )
      {
        ++columnNormalizationIndex;
      }
      double denominatorValue( rotationMatrix( ( columnIndex + 1 ),
                                               referenceField ) );
      squaredSum += nextSquare;
      nextSquare = ( denominatorValue * denominatorValue );
      inverseLengthSquared = ( 1.0
          / ( squaredSum
              * ( 1.0 + ( squaredSum / nextSquare ) ) ) );
      for( size_t rowIndex( 0 );
          rowIndex <= columnIndex;
          ++rowIndex )
      {
        rotationMatrix( rowIndex,
                        columnNormalizationIndex ) *= inverseLengthSquared;
      }
    }
    // So now rotationMatrix is the rotation matrix that takes a vector aligned
    // with the axis of field referenceField to (endNode - startNode).

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "rotationMatrix ="
    << std::endl << rotationMatrix;
    Eigen::VectorXd referenceVector( numberOfFields );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex == referenceField )
      {
        referenceVector( fieldIndex ) = 1.0;
      }
      else
      {
        referenceVector( fieldIndex ) = 0.0;
      }
    }
    std::cout << std::endl
    << "referenceVector ="
    << std::endl << referenceVector;
    std::cout << std::endl
    << "rotationMatrix * referenceVector ="
    << std::endl << ( rotationMatrix * referenceVector );
    std::cout << std::endl
    << "endNode - startNode = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << ( endNode[ fieldIndex ] - startNode[ fieldIndex ] );
    }
    std::cout << " }";
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
