/*
 * NodesOnParallelPlanes.cpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/NodesOnParallelPlanes.hpp"

namespace VevaciousPlusPlus
{

  NodesOnParallelPlanes::NodesOnParallelPlanes(
                                      std::vector< double > const& falseVacuum,
                                       std::vector< double > const& trueVacuum,
                                     size_t const numberOfIntermediateNodes ) :
    NodesOnPlanes( falseVacuum,
                   trueVacuum,
                   numberOfIntermediateNodes )
  {
    // This constructor is just an initialization list.
  }

  NodesOnParallelPlanes::~NodesOnParallelPlanes()
  {
    // This does nothing.
  }


  // This takes nodeParameterization as a vector in the plane with field
  // referenceField = 0 and projects it onto the plane perpendicular to the
  // difference vector between the vacua, and adds that to nodeVector.
  void NodesOnParallelPlanes::AddTransformedNode(
                                       std::vector< double > const& nodeVector,
                                        std::vector< double > const& startNode,
                                          std::vector< double > const& endNode,
                      std::vector< double > const& nodeParameterization ) const
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "NodesOnParallelPlanes::AddTransformedNode( nodeVector = {";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << nodeVector[ fieldIndex ];
    }
    std::cout << " }, startNode = {";
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
    std::cout << " }, endNode = {";
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
    std::cout << " }, nodeParameterization = {";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfParametersPerNode;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << nodeParameterization[ fieldIndex ];
    }
    std::cout << " } ) called.";
    std::cout << std::endl;/**/
    // Shorthand:
    // n = nodeVector
    // d = ( endNode - startNode )
    // p = vector in plane from nodeParameterization
    // n += p - d ( p.d / d.d )
    double dotProduct( 0.0 );
    double differenceSquared( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double const differenceValue( endNode[ fieldIndex ]
                                    - startNode[ fieldIndex ] );
      differenceSquared += ( differenceValue * differenceValue );
      if( fieldIndex < referenceField )
      {
        dotProduct += ( differenceValue * nodeParameterization[ fieldIndex ] );
      }
      else if( fieldIndex > referenceField )
      {
        dotProduct
        += ( differenceValue * nodeParameterization[ fieldIndex - 1 ] );
      }
    }
    double const subtractionFraction( dotProduct / differenceSquared );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      nodeVector[ fieldIndex ] += ( subtractionFraction
                       * ( startNode[ fieldIndex ] - endNode[ fieldIndex ] ) );
      if( fieldIndex < referenceField )
      {
        nodeVector[ fieldIndex ] += nodeParameterization[ fieldIndex ];
      }
      else if( fieldIndex > referenceField )
      {
        nodeVector[ fieldIndex ] += nodeParameterization[ fieldIndex - 1 ];
      }
    }
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
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
    std::cout << std::endl
    << "node in plane = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      if( fieldIndex < referenceField )
      {
        std::cout << nodeParameterization[ fieldIndex ];
      }
      else if( fieldIndex == referenceField )
      {
        std::cout << "0.0";
      }
      else if( fieldIndex > referenceField )
      {
        std::cout << nodeParameterization[ fieldIndex - 1 ];
      }
    }
    std::cout << " }";
    std::cout << std::endl
    << "transformed node = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << nodeVector[ fieldIndex ];
    }
    std::cout << " }";
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
