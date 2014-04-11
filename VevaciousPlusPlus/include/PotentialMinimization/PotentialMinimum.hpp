/*
 * PotentialMinimum.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALMINIMUM_HPP_
#define POTENTIALMINIMUM_HPP_

#include "../StandardIncludes.hpp"
#include "MinuitMinimum.hpp"

namespace VevaciousPlusPlus
{
  class PotentialMinimum : public MinuitMinimum
  {
  public:
    PotentialMinimum( std::vector< double > const& fieldConfiguration,
                      double const potentialDepth );
    PotentialMinimum( MinuitMinimum const& minuitMinimum );
    PotentialMinimum();
    virtual
    ~PotentialMinimum();


    // This returns the sum of the squares of the differences in the field
    // values of this PotentialMinimum with comparisonMinimum.
    double SquareDistanceTo( PotentialMinimum const& comparisonMinimum ) const;

    // This returns the sum of the squares of the field values.
    double LengthSquared() const;

    std::vector< double > const& FieldConfiguration() const
    { return variableValues; }

    double PotentialValue() const{ return functionValue; }
  };




  // This returns the sum of the squares of the differences in the field values
  // of this PotentialMinimum with comparisonMinimum.
  inline double PotentialMinimum::SquareDistanceTo(
                              PotentialMinimum const& comparisonMinimum ) const
  {
    double returnDouble( 0.0 );
    double fieldDifference( 0.0 );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < variableValues.size();
         ++fieldIndex )
    {
      fieldDifference = ( variableValues[ fieldIndex ]
                          - comparisonMinimum.VariableValues()[ fieldIndex ] );
      returnDouble += ( fieldDifference * fieldDifference );
    }
    return returnDouble;
  }

  // This returns the sum of the squares of the field values.
  inline double PotentialMinimum::LengthSquared() const
  {
    double returnDouble( 0.0 );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < variableValues.size();
         ++fieldIndex )
    {
      returnDouble += ( variableValues[ fieldIndex ]
                       * variableValues[ fieldIndex ] );
    }
    return returnDouble;
  }

} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALMINIMUM_HPP_ */
