/*
 * QuadraticSplineThroughNodes.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef QUADRATICSPLINETHROUGHNODES_HPP_
#define QUADRATICSPLINETHROUGHNODES_HPP_

#include "CommonIncludes.hpp"
#include "TunnelPath.hpp"
#include "Eigen/Dense"
#include "QuadraticSplinePathSegment.hpp"

namespace VevaciousPlusPlus
{

  class QuadraticSplineThroughNodes : public TunnelPath
  {
  public:
    QuadraticSplineThroughNodes(
                         std::vector< std::vector< double > > const& pathNodes,
                             std::vector< double > const& pathParameterization,
                                 double const pathTemperature );
    virtual ~QuadraticSplineThroughNodes();


    // This fills fieldConfiguration with the values that the fields
    // should have when the path auxiliary is given by auxiliaryValue.
    void PutOnPathAt( std::vector< double >& fieldConfiguration,
                      double const auxiliaryValue ) const;

    // This returns the dot product with itself of the derivative of the
    // field vector with respect to the path auxiliary evaluated at
    // auxiliaryValue.
    double SlopeSquared( double const auxiliaryValue ) const;

    // This returns the dot product of the first derivative of the field
    // vector with the second derivative, both with respect to the path
    // auxiliary, evaluated at auxiliaryValue.
    double SlopeDotAcceleration( double const auxiliaryValue ) const;

    // This is for debugging.
    std::string AsDebuggingString() const;


  protected:
    size_t const numberOfSegments;
    double const inverseSegmentLength;
    double const segmentLength;
    // The segments are of equal length, so we jump straight to the correct
    // segment for the various functions.
    std::vector< QuadraticSplinePathSegment > pathSegments;


    // This gives the index for which path segment is correct for
    // auxiliaryValue along with the value of the auxiliary value along the
    // segment.
    std::pair< size_t, double >
    SegmentAuxiliary( double const auxiliaryValue ) const;
  };




  // This fills fieldConfiguration with the values that the fields
  // should have when the path auxiliary is given by auxiliaryValue.
  inline void QuadraticSplineThroughNodes::PutOnPathAt(
                                     std::vector< double >& fieldConfiguration,
                                            double const auxiliaryValue ) const
  {
    std::pair< size_t, double > const
    indexAndRemainder( SegmentAuxiliary( auxiliaryValue ) );
    pathSegments[ indexAndRemainder.first ].PutOnSegment( fieldConfiguration,
                                                    indexAndRemainder.second );
  }

  // This returns the dot product with itself of the derivative of the
  // field vector with respect to the path auxiliary evaluated at
  // auxiliaryValue.
  inline double QuadraticSplineThroughNodes::SlopeSquared(
                                            double const auxiliaryValue ) const
  {
    std::pair< size_t, double > const
    indexAndRemainder( SegmentAuxiliary( auxiliaryValue ) );
    return pathSegments[ indexAndRemainder.first ].SlopeSquared(
                                                    indexAndRemainder.second );
  }

  // This returns the dot product of the first derivative of the field
  // vector with the second derivative, both with respect to the path
  // auxiliary, evaluated at auxiliaryValue.
  inline double QuadraticSplineThroughNodes::SlopeDotAcceleration(
                                            double const auxiliaryValue ) const
  {
    std::pair< size_t, double > const
    indexAndRemainder( SegmentAuxiliary( auxiliaryValue ) );
    return pathSegments[ indexAndRemainder.first ].SlopeSquared(
                                                    indexAndRemainder.second );
  }

  // This gives the index for which path segment is correct for
  // auxiliaryValue along with the value of the auxiliary value along the
  // segment.
  inline std::pair< size_t, double >
  QuadraticSplineThroughNodes::SegmentAuxiliary(
                                            double const auxiliaryValue ) const
  {
    size_t segmentIndex( auxiliaryValue * inverseSegmentLength );
    if( auxiliaryValue < 0.0 )
    {
      segmentIndex = 0;
    }
    else if( segmentIndex >= pathSegments.size() )
    {
      segmentIndex = ( pathSegments.size() - 1 );
    }
    return std::make_pair( segmentIndex,
                       ( auxiliaryValue - ( segmentIndex * segmentLength ) ) );
  }

  // This is for debugging.
  inline std::string QuadraticSplineThroughNodes::AsDebuggingString() const
  {
    double const segmentLength( 1.0 / inverseSegmentLength );
    std::stringstream returnStream;
    returnStream << "{ ";
    for( size_t segmentIndex( 0 );
         segmentIndex < pathSegments.size();
         ++segmentIndex )
    {
      if( segmentIndex > 0 )
      {
        returnStream << "," << std::endl;
      }
      returnStream << pathSegments[ segmentIndex ].AsDebuggingString(
                                                segmentIndex * segmentLength );
    }
    returnStream << " }";
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* QUADRATICSPLINETHROUGHNODES_HPP_ */
