/*
 * LinearSplinePathSegment.hpp
 *
 *  Created on: Jul 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LINEARSPLINEPATHSEGMENT_HPP_
#define LINEARSPLINEPATHSEGMENT_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class LinearSplinePathSegment
  {
  public:
    LinearSplinePathSegment( std::vector< double > const& startNode,
                             std::vector< double > const& endNode,
                             double const segmentAuxiliaryLength );
    LinearSplinePathSegment( LinearSplinePathSegment const& copySource );
    LinearSplinePathSegment();
    virtual ~LinearSplinePathSegment();


    // This fills fieldConfiguration with the values that the fields should
    // have when the segment auxiliary is given by segmentAuxiliary.
    void PutOnSegment( std::vector< double >& fieldConfiguration,
                       double const segmentAuxiliary ) const;

    // This returns the sum of the squares of the slopes at segmentAuxiliary.
    double SlopeSquared( double const segmentAuxiliary ) const;

    // This should return the dot product of the first derivative of the field
    // vector with the second derivative, both with respect to the segment
    // auxiliary, evaluated at segmentAuxiliary.
    double SlopeDotAcceleration( double const segmentAuxiliary ) const
    { return 0.0; }

    // This returns the slope at the end of the segment, which in this case is
    // the same as the slope everywhere in the segment.
    std::vector< double > EndSlope() const{ return fieldLinears; }

    double SegmentLength() const{ return segmentAuxiliaryLength; }

    // This is for debugging.
    std::string AsDebuggingString( double const segmentStart ) const;


  protected:
    size_t numberOfFields;
    std::vector< double > fieldConstants;
    std::vector< double > fieldLinears;
    double segmentAuxiliaryLength;
  };




  // This fills fieldConfiguration with the values that the fields should
  // have when the segment auxiliary is given by segmentAuxiliary.
  inline void LinearSplinePathSegment::PutOnSegment(
                                     std::vector< double >& fieldConfiguration,
                                          double const segmentAuxiliary ) const
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      fieldConfiguration[ fieldIndex ]
      = ( fieldConstants[ fieldIndex ]
          + ( segmentAuxiliary * fieldLinears[ fieldIndex ] ) );
    }
  }

  // This returns the sum of the squares of the slopes at segmentAuxiliary.
  inline double LinearSplinePathSegment::SlopeSquared(
                                          double const segmentAuxiliary ) const
  {
    double returnValue( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double const slopeValue( fieldLinears[ fieldIndex ] );
      returnValue += ( slopeValue * slopeValue );
    }
    return returnValue;
  }

  // This is for debugging.
  inline std::string LinearSplinePathSegment::AsDebuggingString(
                                              double const segmentStart ) const
  {
    std::stringstream returnStream;
    returnStream << "{ ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        returnStream << "," << std::endl << " ";
      }
      returnStream << "( UnitStep[ x - " << segmentStart << " ] * ( "
      << fieldConstants[ fieldIndex ] << " + ( " << fieldLinears[ fieldIndex ]
      << " * ( x - " << segmentStart << " ) ) ) * UnitStep[ "
      << ( segmentStart + segmentAuxiliaryLength ) << " - x ] )";
    }
    returnStream << " }";
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* LINEARSPLINEPATHSEGMENT_HPP_ */
