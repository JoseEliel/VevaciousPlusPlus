/*
 * PathFieldsAndPotential.hpp
 *
 *  Created on: Jun 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PATHFIELDSANDPOTENTIAL_HPP_
#define PATHFIELDSANDPOTENTIAL_HPP_

#include "../CommonIncludes.hpp"
#include "../PotentialEvaluation/SimplePolynomial.hpp"
#include "Eigen/Dense"
#include "SplinePotential.hpp"

namespace VevaciousPlusPlus
{

  class PathFieldsAndPotential
  {
  public:
    PathFieldsAndPotential( Eigen::MatrixXd const& pathCoefficients,
                         std::vector< double > const& falseVacuumConfiguration,
                            double const falseVacuumDepth,
                            double const trueVacuumDepth,
                            double const givenTemperature );
    PathFieldsAndPotential( PathFieldsAndPotential const& copySource );
    virtual
    ~PathFieldsAndPotential();


    double TruePotential() const{ return trueVacuumDepth; }

    double FalsePotential() const{ return falseVacuumDepth; }

    double PotentialApproximation( double const auxiliaryValue ) const
    { return potentialApproximation( auxiliaryValue ); }

    SplinePotential& PotentialApproximation(){ return potentialApproximation; }
    SplinePotential const& PotentialApproximation() const
    { return potentialApproximation; }

    std::vector< SimplePolynomial > const& FieldPath() const
    { return fieldPath; }

    // This evaluates fieldPath at auxiliaryValue, putting it in
    // fieldConfiguration, and returns fieldConfiguration.
    std::vector< double > const&
    FieldConfiguration( double const auxiliaryValue );

    double FieldDerivative( size_t const fieldIndex,
                            double const auxiliaryValue ) const
    { return pathTangent[ fieldIndex ]( auxiliaryValue ); }

    std::vector< SimplePolynomial > const& FieldDerivatives() const
    { return pathTangent; }

    // This returns the sum of the squares of the field derivatives evaluated
    // at auxiliaryValue.
    double FieldDerivativesSquared( double const auxiliaryValue ) const;

    bool NonZeroTemperature() const{ return nonZeroTemperature; }

    double GivenTemperature() const{ return givenTemperature; }

    // This returns the dimensionality of the radial integral: 2.0 for non-zero
    // temperature, 3.0 for zero temperature.
    double DampingFactor() const;

    // This is for debugging.
    std::string AsDebuggingString() const;

    // This is for debugging.
    std::string FieldsString( double const auxiliaryValue ) const;


  protected:
    SplinePotential potentialApproximation;
    size_t const numberOfFields;
    std::vector< SimplePolynomial > fieldPath;
    std::vector< SimplePolynomial > pathTangent;
    std::vector< double > fieldConfiguration;
    double falseVacuumDepth;
    double trueVacuumDepth;
    bool nonZeroTemperature;
    double givenTemperature;
  };




  // This evaluates fieldPath at auxiliaryValue, putting it in
  // fieldConfiguration, and returns fieldConfiguration.
  inline std::vector< double > const&
  PathFieldsAndPotential::FieldConfiguration( double const auxiliaryValue )
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < fieldPath.size();
         ++fieldIndex )
    {
      fieldConfiguration[ fieldIndex ]
      = fieldPath[ fieldIndex ]( auxiliaryValue );
    }
    return fieldConfiguration;
  }

  // This returns the sum of the squares of the field derivatives evaluated
  // at auxiliaryValue.
  inline double PathFieldsAndPotential::FieldDerivativesSquared(
                                            double const auxiliaryValue ) const
  {
    double returnValue( 0.0 );
    double derivativeValue( NAN );
    for( size_t fieldIndex( 0 );
         fieldIndex < pathTangent.size();
         ++fieldIndex )
    {
      derivativeValue = pathTangent[ fieldIndex ]( auxiliaryValue );
      returnValue += ( derivativeValue * derivativeValue );
    }
    return returnValue;
  }

  // This returns the dimensionality of the radial integral: 2.0 for non-zero
  // temperature, 3.0 for zero temperature.
  inline double PathFieldsAndPotential::DampingFactor() const
  {
    if( nonZeroTemperature )
    {
      return 2.0;
    }
    else
    {
      return 3.0;
    }
    // I could have written
    // return ( nonZeroTemperature ? 2.0 : 3.0 );
    // but I feel that it's nice to be verbose.
  }

  // This is for debugging.
  inline std::string
  PathFieldsAndPotential::FieldsString( double const auxiliaryValue ) const
  {
    std::stringstream returnStream;
    for( size_t fieldIndex( 0 );
         fieldIndex < fieldPath.size();
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        returnStream << ", ";
      }
      returnStream << "f[" << fieldIndex << "] = "
      << fieldPath[ fieldIndex ]( auxiliaryValue );
    }
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* PATHFIELDSANDPOTENTIAL_HPP_ */
