/*
 * QuadraticSplinePathSegment.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/QuadraticSplinePathSegment.hpp"

namespace VevaciousPlusPlus
{

  QuadraticSplinePathSegment::QuadraticSplinePathSegment(
                                        std::vector< double > const& startNode,
                                       std::vector< double > const& startSlope,
                                          std::vector< double > const& endNode,
                                        double const segmentAuxiliaryLength ) :
    numberOfFields( startNode.size() ),
    fieldConstants( startNode ),
    fieldLinears( startSlope ),
    fieldQuadratics( numberOfFields ),
    segmentAuxiliaryLength( segmentAuxiliaryLength )
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      fieldQuadratics[ fieldIndex ]
      = ( ( endNode[ fieldIndex ] - startNode[ fieldIndex ]
            - ( segmentAuxiliaryLength * startSlope[ fieldIndex ] ) )
          / ( segmentAuxiliaryLength * segmentAuxiliaryLength ) );
    }
  }

  // This constructor makes just a straight line.
  QuadraticSplinePathSegment::QuadraticSplinePathSegment(
                                        std::vector< double > const& startNode,
                                          std::vector< double > const& endNode,
                                        double const segmentAuxiliaryLength ) :
    numberOfFields( startNode.size() ),
    fieldConstants( startNode ),
    fieldLinears( numberOfFields ),
    fieldQuadratics( numberOfFields ),
    segmentAuxiliaryLength( segmentAuxiliaryLength )
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      fieldLinears[ fieldIndex ]
      = ( ( endNode[ fieldIndex ] - startNode[ fieldIndex ] )
          / segmentAuxiliaryLength );
      fieldQuadratics[ fieldIndex ] = 0.0;
    }
  }

  QuadraticSplinePathSegment::QuadraticSplinePathSegment(
                               QuadraticSplinePathSegment const& copySource ) :
    numberOfFields( copySource.numberOfFields ),
    fieldConstants( copySource.fieldConstants ),
    fieldLinears( copySource.fieldLinears ),
    fieldQuadratics( copySource.fieldQuadratics ),
    segmentAuxiliaryLength( copySource.segmentAuxiliaryLength )
  {
    // This constructor is just an initialization list.
  }

  QuadraticSplinePathSegment::QuadraticSplinePathSegment() :
    numberOfFields( 0 ),
    fieldConstants(),
    fieldLinears(),
    fieldQuadratics(),
    segmentAuxiliaryLength( NAN )
  {
    // This constructor is just an initialization list.
  }

  QuadraticSplinePathSegment::~QuadraticSplinePathSegment()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
