/*
 * LinearSplineThroughNodes.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LINEARSPLINETHROUGHNODES_HPP_
#define LINEARSPLINETHROUGHNODES_HPP_

#include "CommonIncludes.hpp"
#include "TunnelPath.hpp"
#include "LinearSplinePathSegment.hpp"

namespace VevaciousPlusPlus
{

  class LinearSplineThroughNodes : public TunnelPath
  {
  public:
    LinearSplineThroughNodes(
                         std::vector< std::vector< double > > const& pathNodes,
                             std::vector< double > const& pathParameterization,
                              double const pathTemperature );
    virtual ~LinearSplineThroughNodes();


    // This fills fieldConfiguration with the values that the fields
    // should have when the path auxiliary is given by auxiliaryValue.
    void PutOnPathAt( std::vector< double >& fieldConfiguration,
                      double const auxiliaryValue ) const;

    // This returns the dot product with itself of the derivative of the
    // field vector with respect to the path auxiliary evaluated at
    // auxiliaryValue. This is constant by construction (constant "speed" along
    // straight segments with infinitesimal circle segments between straight
    // segments, where the "acceleration" is perpendicular to the "velocity".)
    double SlopeSquared( double const auxiliaryValue ) const
    { return slopeSquared; }

    // This returns the dot product of the first derivative of the field
    // vector with the second derivative, both with respect to the path
    // auxiliary, evaluated at auxiliaryValue. This is 0.0 by construction
    // (constant "speed" along straight segments with infinitesimal circle
    // segments between straight segments, where the "acceleration" is
    // perpendicular to the "velocity".)
    double SlopeDotAcceleration( double const auxiliaryValue ) const
    { return 0.0; }

    // This is for debugging.
    std::string AsDebuggingString() const;


  protected:
    std::vector< LinearSplinePathSegment > pathSegments;
    double slopeSquared;


    // This gives the index for which path segment is correct for
    // auxiliaryValue along with the value of the auxiliary value along the
    // segment.
    std::pair< size_t, double >
    SegmentAuxiliary( double const auxiliaryValue ) const;
  };




  // This fills fieldConfiguration with the values that the fields
  // should have when the path auxiliary is given by auxiliaryValue.
  inline void LinearSplineThroughNodes::PutOnPathAt(
                                     std::vector< double >& fieldConfiguration,
                                            double const auxiliaryValue ) const
  {
    size_t segmentIndex( 0 );
    double segmentAuxiliary( auxiliaryValue );
    while( ( segmentIndex < ( pathSegments.size() - 1 ) )
           &&
           ( segmentAuxiliary
             > pathSegments[ segmentIndex ].SegmentLength() ) )
    {
      segmentAuxiliary -= pathSegments[ segmentIndex ].SegmentLength();
      ++segmentIndex;
    }
    // At this point, either we found the segment which contains the
    // auxiliary point, or segmentIndex = ( pathSegments.size() - 1 ), so the
    // last segment. Auxiliary values outside the range get put on the
    // extensions of the first or last segment appropriately.
    pathSegments[ segmentIndex ].PutOnSegment( fieldConfiguration,
                                               segmentAuxiliary );
  }

  // This gives the index for which path segment is correct for
  // auxiliaryValue along with the value of the auxiliary value along the
  // segment.
  inline std::pair< size_t, double >
  LinearSplineThroughNodes::SegmentAuxiliary(
                                            double const auxiliaryValue ) const
  {
    if( auxiliaryValue < 0.0 )
    {
      return std::make_pair( 0,
                             0.0 );
    }
    size_t segmentIndex( 0 );
    double auxiliaryRemainder( auxiliaryValue );
    while( segmentIndex < pathSegments.size() )
    {
      if( auxiliaryRemainder < pathSegments[ segmentIndex ].SegmentLength() )
      {
        return std::make_pair( segmentIndex,
                               auxiliaryRemainder );
      }
      auxiliaryRemainder -= pathSegments[ segmentIndex ].SegmentLength();
      ++segmentIndex;
    }
    return std::make_pair( ( pathSegments.size() - 1 ),
                           pathSegments.back().SegmentLength() );
  }

  // This is for debugging.
  inline std::string LinearSplineThroughNodes::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream << "{ ";
    double cumulativeAuxiliary( 0.0 );
    for( size_t segmentIndex( 0 );
         segmentIndex < pathSegments.size();
         ++segmentIndex )
    {
      if( segmentIndex > 0 )
      {
        returnStream << "," << std::endl;
      }
      returnStream
      << pathSegments[ segmentIndex ].AsDebuggingString( cumulativeAuxiliary );
      cumulativeAuxiliary += pathSegments[ segmentIndex ].SegmentLength();
    }
    returnStream << " }";
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* LINEARSPLINETHROUGHNODES_HPP_ */
