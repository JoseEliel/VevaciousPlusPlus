/*
 * LinearSplineThroughNodes.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/LinearSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  LinearSplineThroughNodes::LinearSplineThroughNodes(
                         std::vector< std::vector< double > > const& pathNodes,
                             std::vector< double > const& pathParameterization,
                                               double const pathTemperature ) :
    TunnelPath( pathNodes.front().size(),
                pathParameterization,
                pathTemperature ),
    pathSegments( pathNodes.size() - 1 ),
    slopeSquared( -1.0 )
  {
    std::vector< double > segmentLengths( pathSegments.size() );
    double totalLength( 0.0 );
    for( size_t segmentIndex( 0 );
         segmentIndex < pathSegments.size();
         ++segmentIndex )
    {
      double squaredLength( 0.0 );
      std::vector< double > const& startNode( pathNodes[ segmentIndex ] );
      std::vector< double > const& endNode( pathNodes[ segmentIndex + 1 ] );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        double const
        fieldDifference( endNode[ fieldIndex ] - startNode[ fieldIndex ] );
        squaredLength += ( fieldDifference * fieldDifference );
      }
      segmentLengths[ segmentIndex ] = sqrt( squaredLength );
      totalLength += segmentLengths[ segmentIndex ];
    }

    // The constant velocity is totalLength divided by the range of the path
    // auxiliary, which goes from 0 to 1, so the dot product of the vector of
    // the derivatives of the fields with respect to the path auxiliary with
    // itself is just ( totalLength / 1.0 )^2.
    slopeSquared = ( totalLength * totalLength );
    double const inverseTotalLength( 1.0 / totalLength );
    for( size_t segmentIndex( 0 );
         segmentIndex < pathSegments.size();
         ++segmentIndex )
    {
      pathSegments[ segmentIndex ]
      = LinearSplinePathSegment( pathNodes[ segmentIndex ],
                                 pathNodes[ segmentIndex + 1 ],
                     ( segmentLengths[ segmentIndex ] * inverseTotalLength ) );
    }
  }

  LinearSplineThroughNodes::~LinearSplineThroughNodes()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
