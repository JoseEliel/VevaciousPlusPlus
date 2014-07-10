/*
 * QuadraticSplineThroughNodes.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/QuadraticSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  QuadraticSplineThroughNodes::QuadraticSplineThroughNodes(
                         std::vector< std::vector< double > > const& pathNodes,
                                               double const pathTemperature ) :
    TunnelPath( pathNodes.front().size(),
                pathTemperature ),
    numberOfSegments( pathNodes.size() - 1 ),
    inverseSegmentLength( 1.0 / (double)numberOfSegments ),
    pathSegments( numberOfSegments )
  {
    double const segmentLength( (double)numberOfSegments );
    pathSegments.front() = QuadraticSplinePathSegment( pathNodes.front(),
                                                       pathNodes[ 1 ],
                                                       segmentLength );
    for( size_t segmentIndex( 1 );
         segmentIndex < numberOfSegments;
         ++segmentIndex )
    {
      pathSegments[ segmentIndex ]
      = QuadraticSplinePathSegment( pathNodes[ segmentIndex ],
                                   pathSegments[ segmentIndex - 1 ].EndSlope(),
                                    pathNodes[ segmentIndex + 1 ],
                                    segmentLength );
    }
  }

  QuadraticSplineThroughNodes::~QuadraticSplineThroughNodes()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
