/*
 * QuadraticSplineThroughNodes.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/QuadraticSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  QuadraticSplineThroughNodes::QuadraticSplineThroughNodes(
                         std::vector< std::vector< double > > const& pathNodes,
                             std::vector< double > const& pathParameterization,
                                               double const pathTemperature ) :
    TunnelPath( pathNodes.front().size(),
                pathParameterization,
                pathTemperature ),
    numberOfSegments( pathNodes.size() - 1 ),
    inverseSegmentLength( (double)numberOfSegments ),
    pathSegments( numberOfSegments )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "QuadraticSplineThroughNodes::QuadraticSplineThroughNodes( pathNodes ="
    << " { ";
    for( size_t nodeIndex( 0 );
         nodeIndex < pathNodes.size();
         ++nodeIndex )
    {
      if( nodeIndex > 0 )
      {
        std::cout << " }," << std::endl;
      }
      std::cout << "{ ";
      for( size_t fieldIndex( 0 );
           fieldIndex < pathNodes[ nodeIndex ].size();
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << pathNodes[ nodeIndex ][ fieldIndex ];
      }
    }
    std::cout << " } }, pathParameterization = { " << std::endl;
    for( std::vector< double >::const_iterator
         pathParameter( pathParameterization.begin() );
         pathParameter < pathParameterization.end();
         ++pathParameter )
    {
      if( pathParameter > pathParameterization.begin() )
      {
        std::cout << ", ";
      }
      std::cout << *pathParameter;
    }
    std::cout << " }, pathTemperature = " << pathTemperature
    << " ) called.";
    std::cout << std::endl;/**/

    double const segmentLength( 1.0 / (double)numberOfSegments );
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

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "this->AsDebuggingString() = " << this->AsDebuggingString();
    std::cout << std::endl;/**/
  }

  QuadraticSplineThroughNodes::~QuadraticSplineThroughNodes()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
