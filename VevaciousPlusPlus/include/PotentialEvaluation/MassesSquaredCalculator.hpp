/*
 * MassesSquaredCalculator.hpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MASSESSQUAREDCALCULATOR_HPP_
#define MASSESSQUAREDCALCULATOR_HPP_

#include <map>
#include <string>
#include <vector>
#include "LHPC/Utilities/ParsingUtilities.hpp"
#include <stdexcept>

namespace VevaciousPlusPlus
{

  class MassesSquaredCalculator
  {
  public:
    enum SpinType
    {
      scalarBoson,
      weylFermion,
      gaugeBoson,
      notSet
    };

    MassesSquaredCalculator(
                    std::map< std::string, std::string > const& attributeMap );

    MassesSquaredCalculator( MassesSquaredCalculator const& copySource ) :
      multiplicityFactor( copySource.multiplicityFactor ),
      spinType( copySource.spinType ) {}

    MassesSquaredCalculator() :
      multiplicityFactor( 0.0 ),
      spinType( notSet ) {}

    virtual ~MassesSquaredCalculator() {}


    // This should return the masses-squared using the values for the
    // Lagrangian parameters found in parameterValues and the values for the
    // fields found in fieldConfiguration.
    virtual std::vector< double >
    MassesSquared( std::vector< double > const& parameterValues,
                   std::vector< double > const& fieldConfiguration ) const = 0;

    // This should return the masses-squared using the values for the
    // Lagrangian parameters from the last call of UpdateForFixedScale and the
    // values for the fields found in fieldConfiguration.
    virtual std::vector< double >
    MassesSquared( std::vector< double > const& fieldConfiguration ) const = 0;

    // This returns the number of identical copies of this mass-squared matrix
    // that the model has.
    double MultiplicityFactor() const{ return multiplicityFactor; }

    SpinType GetSpinType() const{ return spinType; }

    // This should update all objects which contribute to the masses with the
    // values for the Lagrangian parameters given in parameterValues.
    virtual void
    UpdateForFixedScale( std::vector< double > const& parameterValues ) = 0;


  protected:
    double multiplicityFactor;
    SpinType spinType;


    // This sets the spin type based on the string found for the attribute. It
    // throws an exception if it is not a valid spin type.
    void SetSpinType( std::string const& attributeValue );
  };





  inline MassesSquaredCalculator::MassesSquaredCalculator(
                   std::map< std::string, std::string > const& attributeMap ) :
    multiplicityFactor( 1.0 ),
    spinType( notSet )
  {
    std::map< std::string, std::string >::const_iterator
    attributeFinder( attributeMap.find( "MultiplicityFactor" ) );
    if( attributeFinder != attributeMap.end() )
    {
      multiplicityFactor
      = LHPC::ParsingUtilities::StringToDouble( attributeFinder->second );
    }
    attributeFinder = attributeMap.find( "SpinType" );
    if( attributeFinder == attributeMap.end() )
    {
      throw std::runtime_error(
                    "Mass(-squared) matrix needs an attribute \"SpinType\"!" );
    }
    SetSpinType( attributeFinder->second );
  }

  // This sets the spin type based on the string found for the attribute. It
  // throws an exception if it is not a valid spin type.
  inline void
  MassesSquaredCalculator::SetSpinType( std::string const& attributeValue )
  {
    if( attributeValue == "ScalarBoson" )
    {
      spinType = scalarBoson;
    }
    else if( attributeValue == "WeylFermion" )
    {
      spinType = weylFermion;
    }
    else if(attributeValue == "GaugeBoson" )
    {
      spinType = gaugeBoson;
    }
    else
    {
      std::stringstream errorBuilder;
      errorBuilder << "Invalid value for SpinType (\"" << attributeValue
      << "\"). SpinType must be \"ScalarBoson\", \"WeylFermion\", or"
      << " \"GaugeBoson\":";
      throw std::runtime_error( errorBuilder.str() );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* MASSESSQUAREDCALCULATOR_HPP_ */
