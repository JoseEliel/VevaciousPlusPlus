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
    reflectionMatrices()
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

    // We ensure that pathNodes and reflectionMatrices have the correct size.
    pathNodes.resize( this->numberOfIntermediateNodes + 2 );
    sideNodeIndices.resize( pathNodes.size() );
    reflectionMatrices.resize( pathNodes.size(),
                             Eigen::MatrixXd( numberOfFields,
                                              numberOfFields ) );
    // There are only rotation matrices for the variable nodes.

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "NodesOnBisectingPlanes::NodesOnBisectingPlanes( numberOfFields = "
    << numberOfFields << ", numberOfIntermediateNodes = "
    << numberOfIntermediateNodes << " ) just made rotationMatrices: {"
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

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "NodesOnBisectingPlanes::NodesOnBisectingPlanes( numberOfFields = "
    << numberOfFields << ", numberOfIntermediateNodes = "
    << numberOfIntermediateNodes << " ) finished."
    << std::endl;
    std::cout << "adjustmentOrder = { ";
    for( size_t orderIndex( 0 );
         orderIndex < adjustmentOrder.size();
         ++orderIndex )
    {
      if( orderIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << adjustmentOrder[ orderIndex ];
    }
    std::cout << "}" << std::endl;
    std::cout << "sideNodeIndices = { ";
    for( size_t sideIndex( 0 );
         sideIndex < sideNodeIndices.size();
         ++sideIndex )
    {
      if( sideIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << "[ " << sideNodeIndices[ sideIndex ].first << ", "
      << sideNodeIndices[ sideIndex ].second << " ]";
    }
    std::cout << "}" << std::endl;
    std::cout << "pathNodes = { { ";
    for( size_t nodeIndex( 0 );
         nodeIndex < pathNodes.size();
         ++nodeIndex )
    {
      if( nodeIndex > 0 )
      {
        std::cout << " }," << std::endl << "{ ";
      }
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << pathNodes[ nodeIndex ][ fieldIndex ];
      }
    }
    std::cout << " } }" << std::endl;
    std::cout << std::endl;/**/
  }

  NodesOnBisectingPlanes::~NodesOnBisectingPlanes()
  {
    // This does nothing.
  }


  // This sets reflectionMatrices[ nodeIndex ] to be a matrix that rotates the
  // vector difference from FalseSideNode( nodeIndex, pathNodes ) ) to
  // TrueSideNode( nodeIndex, pathNodes ) to lie along the axis of
  // referenceField.
  void NodesOnBisectingPlanes::UpdateRotationMatrix( size_t const nodeIndex )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "NodesOnBisectingPlanes::UpdateRotationMatrix( nodeIndex = "
    << nodeIndex << " ) called.";
    std::cout << std::endl;/**/

    std::vector< double > vacuumDifference;
    SetAsVectorDifference( vacuumDifference,
                           FalseSideNode( nodeIndex,
                                          pathNodes ),
                           TrueSideNode( nodeIndex,
                                         pathNodes ) );
    SetAsHouseholderReflectionFromAxisToVector(
                                               reflectionMatrices[ nodeIndex ],
                                                0,
                                                vacuumDifference );
    // For convenience, the 0th field is taken to be the 0.0 in the
    // parameterization. Taking anything else would complicate the rest of the
    // functions.


    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl;
    std::vector< double > const& startNode( FalseSideNode( nodeIndex,
                                                           pathNodes ) );
    std::vector< double > const& endNode( TrueSideNode( nodeIndex,
                                                        pathNodes ) );
    Eigen::MatrixXd& reflectionMatrix( reflectionMatrices[ nodeIndex ] );
    std::cout << "reflectionMatrix =" << std::endl << reflectionMatrix;
    Eigen::VectorXd referenceVector( Eigen::VectorXd::Zero( numberOfFields ) );
    referenceVector( 0 ) = 1.0;
    std::cout << std::endl
    << "referenceVector ="
    << std::endl << referenceVector;
    std::cout << std::endl
    << "reflectionMatrix * referenceVector ="
    << std::endl << ( reflectionMatrix * referenceVector );
    std::cout << std::endl
    << "startNode = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << startNode[ fieldIndex ];
    }
    std::cout << " }";
    std::cout << std::endl
    << "endNode = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << endNode[ fieldIndex ];
    }
    std::cout << " }";
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
