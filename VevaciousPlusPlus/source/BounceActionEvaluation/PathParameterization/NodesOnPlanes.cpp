/*
 * NodesOnPlanes.cpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/NodesOnPlanes.hpp"

namespace VevaciousPlusPlus
{

  NodesOnPlanes::NodesOnPlanes( size_t const numberOfFields,
                                size_t const numberOfIntermediateNodes ) :
    NodesFromParameterization( numberOfFields,
                               numberOfIntermediateNodes ),
    referenceField( 0 ),
    numberOfParametersPerNode( numberOfFields - 1 )
  {
    // This constructor is just an initialization list.
  }

  NodesOnPlanes::~NodesOnPlanes()
  {
    // This does nothing.
  }


  // This resets the NodesFromParameterization so that it will produce
  // TunnelPath*s that parameterize the path between the given vacua.
  void NodesOnPlanes::SetVacua( PotentialMinimum const& falseVacuum,
                                PotentialMinimum const& trueVacuum )
  {
    // We choose referenceField to be the index of the field with largest
    // magnitude of difference between the vacua.
    double largestDifference( trueVacuum.VariableValues().front()
                              - falseVacuum.VariableValues().front() );
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
      = ( trueVacuum.VariableValues()[ fieldIndex ]
          - falseVacuum.VariableValues()[ fieldIndex ] );
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
    pathNodes.front() = falseVacuum.FieldConfiguration();
    pathNodes.back() = trueVacuum.FieldConfiguration();
    SetInitialParameterizationAndStepSizes( zeroParameterization,
                                            initialStepSizes );
  }

} /* namespace VevaciousPlusPlus */
