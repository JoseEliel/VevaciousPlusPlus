/*
 * NodesOnPlanes.cpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/NodesOnPlanes.hpp"

namespace VevaciousPlusPlus
{

  NodesOnPlanes::NodesOnPlanes( std::vector< double > const& falseVacuum,
                                std::vector< double > const& trueVacuum,
                                size_t const numberOfIntermediateNodes ) :
    NodesFromParameterization( trueVacuum.size(),
                               numberOfIntermediateNodes ),
    referenceField( 0 ),
    numberOfParametersPerNode( numberOfFields - 1 )
  {
    // We choose referenceField to be the index of the field with largest
    // magnitude of difference between the vacua.
    double largestDifference( trueVacuum.front() - falseVacuum.front() );
    if( largestDifference < 0.0 )
    {
      largestDifference = -largestDifference;
    }
    double currentDifference( largestDifference );
    for( size_t fieldIndex( 1 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentDifference
      = ( trueVacuum[ fieldIndex ] - falseVacuum[ fieldIndex ] );
      if( currentDifference < 0.0 )
      {
        currentDifference = -currentDifference;
      }
      if( currentDifference > largestDifference )
      {
        referenceField = fieldIndex;
        largestDifference = currentDifference;
      }
    }
  }

  NodesOnPlanes::~NodesOnPlanes()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
