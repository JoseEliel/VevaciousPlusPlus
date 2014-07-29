/*
 * QuadraticSplinePathSegment.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef QUADRATICSPLINEPATHSEGMENT_HPP_
#define QUADRATICSPLINEPATHSEGMENT_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class QuadraticSplinePathSegment
  {
  public:
    QuadraticSplinePathSegment( std::vector< double > const& startNode,
                                std::vector< double > const& startSlope,
                                std::vector< double > const& endNode,
                                double const segmentAuxiliaryLength );
    QuadraticSplinePathSegment( std::vector< double > const& startNode,
                                std::vector< double > const& endNode,
                                double const segmentAuxiliaryLength );
    QuadraticSplinePathSegment( QuadraticSplinePathSegment const& copySource );
    QuadraticSplinePathSegment();
    virtual ~QuadraticSplinePathSegment();


    // This fills fieldConfiguration with the values that the fields should
    // have when the segment auxiliary is given by segmentAuxiliary.
    void PutOnSegment( std::vector< double >& fieldConfiguration,
                       double const segmentAuxiliary ) const;

    // This returns the sum of the squares of the slopes at segmentAuxiliary.
    double SlopeSquared( double const segmentAuxiliary ) const;

    // This should return the dot product of the first derivative of the field
    // vector with the second derivative, both with respect to the segment
    // auxiliary, evaluated at segmentAuxiliary.
    double SlopeDotAcceleration( double const segmentAuxiliary ) const;

    // This returns the slope at the end of the segment.
    std::vector< double > EndSlope() const;

    double SegmentLength() const{ return segmentAuxiliaryLength; }

    // This is for debugging.
    std::string AsDebuggingString( double const segmentStart ) const;


  protected:
    size_t numberOfFields;
    std::vector< double > fieldConstants;
    std::vector< double > fieldLinears;
    std::vector< double > fieldQuadratics;
    double segmentAuxiliaryLength;
  };




  // This fills fieldConfiguration with the values that the fields should
  // have when the segment auxiliary is given by segmentAuxiliary.
  inline void QuadraticSplinePathSegment::PutOnSegment(
                                     std::vector< double >& fieldConfiguration,
                                          double const segmentAuxiliary ) const
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      fieldConfiguration[ fieldIndex ]
      = ( fieldConstants[ fieldIndex ]
          + ( segmentAuxiliary * ( fieldLinears[ fieldIndex ]
                                   + ( segmentAuxiliary
                                       * fieldQuadratics[ fieldIndex ] ) ) ) );
    }
  }

  // This returns the sum of the squares of the slopes at segmentAuxiliary.
  inline double QuadraticSplinePathSegment::SlopeSquared(
                                          double const segmentAuxiliary ) const
  {
    double returnValue( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double const slopeValue( fieldLinears[ fieldIndex ]
                + ( 2.0 * segmentAuxiliary * fieldQuadratics[ fieldIndex ] ) );
      returnValue += ( slopeValue * slopeValue );
    }
    return returnValue;
  }

  // This should return the dot product of the first derivative of the field
  // vector with the second derivative, both with respect to the segment
  // auxiliary, evaluated at segmentAuxiliary.
  inline double QuadraticSplinePathSegment::SlopeDotAcceleration(
                                          double const segmentAuxiliary ) const
  {
    double returnValue( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double const slopeValue( fieldLinears[ fieldIndex ]
                + ( 2.0 * segmentAuxiliary * fieldQuadratics[ fieldIndex ] ) );
      returnValue += ( slopeValue * 2.0 * fieldQuadratics[ fieldIndex ] );
    }
    return returnValue;
  }

  // This returns the slope at the end of the segment.
  inline std::vector< double > QuadraticSplinePathSegment::EndSlope() const
  {
    std::vector< double > returnVector( fieldLinears );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      returnVector[ fieldIndex ]
      += ( 2.0 * segmentAuxiliaryLength * fieldQuadratics[ fieldIndex ] );
    }
    return returnVector;
  }

  // This is for debugging.
  inline std::string QuadraticSplinePathSegment::AsDebuggingString(
                                              double const segmentStart ) const
  {
    std::stringstream returnStream;
    returnStream << "( UnitStep[ x - " << segmentStart << " ] * ( ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        returnStream << "," << std::endl;
      }
      returnStream << fieldConstants[ fieldIndex ] << " + ( "
      << fieldLinears[ fieldIndex ] << " * ( x - " << segmentStart
      << " ) ) + ( " << fieldQuadratics[ fieldIndex ] << " * ( x - "
      << segmentStart << " )^2 ) ) * UnitStep[ "
      << ( segmentStart + segmentAuxiliaryLength ) << " - x ] )";
    }
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* QUADRATICSPLINEPATHSEGMENT_HPP_ */
