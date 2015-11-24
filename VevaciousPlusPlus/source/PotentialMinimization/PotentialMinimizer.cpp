/*
 * PotentialMinimizer.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/PotentialMinimizer.hpp"

namespace VevaciousPlusPlus
{

  PotentialMinimizer::PotentialMinimizer(
                                       PotentialFunction& potentialFunction ) :
    potentialFunction( potentialFunction ),
    foundMinima(),
    dsbVacuum(),
    panicVacua(),
    panicVacuum()
  {
    // This constructor is just an initialization list.
  }

  PotentialMinimizer::~PotentialMinimizer()
  {
    // This does nothing.
  }


  // This returns false if the vector of the absolute values of the fields of
  // comparisonMinimum lies within a hypercube of side thresholdDistance
  // centered on dsbVacuum.
  bool PotentialMinimizer::IsNotPhaseRotationOfDsbVacuum(
                                     PotentialMinimum const& comparisonMinimum,
                                         double const thresholdDistance ) const
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < dsbVacuum.FieldConfiguration().size();
         ++fieldIndex )
    {
      double dsbField( dsbVacuum.FieldConfiguration()[ fieldIndex ] );
      if( dsbField < 0.0 )
      {
        dsbField = -dsbField;
      }
      double
      comparisonField( comparisonMinimum.FieldConfiguration()[ fieldIndex ] );
      if( comparisonField < 0.0 )
      {
        comparisonField = -comparisonField;
      }
      double fieldDifference( dsbField - comparisonField );
      if( fieldDifference < 0.0 )
      {
        fieldDifference = -fieldDifference;
      }
      if( fieldDifference > thresholdDistance )
      {
        return true;
      }
    }
    return false;
  }

} /* namespace VevaciousPlusPlus */
