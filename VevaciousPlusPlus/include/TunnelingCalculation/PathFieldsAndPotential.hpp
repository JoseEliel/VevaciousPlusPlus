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


    // This sets potentialApproximation to be a SimplePolynomial based on
    // potentialCoefficients with a leading power of 2, also adding a final
    // coefficient (for a term with power 1 higher than given by
    // potentialCoefficients + 2 for the leading power, so its size + 2 (as it
    // gives coefficients for powers 2 to its size + 1)) which makes the
    // potential have a minimum when the auxiliary value is 1.0, which should
    // have been taken into account already in creating potentialCoefficients.
    void SetPotential( Eigen::VectorXd const& potentialCoefficients );

    double TruePotential() const{ return trueVacuumDepth; }

    double FalsePotential() const{ return falseVacuumDepth; }

    double PotentialApproximation( double const auxiliaryValue ) const
    { return potentialApproximation( auxiliaryValue ); }

    SimplePolynomial const& PotentialApproximation() const
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


  protected:
    SimplePolynomial potentialApproximation;
    size_t const numberOfFields;
    std::vector< SimplePolynomial > fieldPath;
    std::vector< SimplePolynomial > pathTangent;
    std::vector< double > fieldConfiguration;
    double falseVacuumDepth;
    double trueVacuumDepth;
    bool nonZeroTemperature;
    double givenTemperature;
  };




  // This sets potentialApproximation to be a SimplePolynomial based on
  // potentialCoefficients with a leading power of 2, also adding a final
  // coefficient (for a term with power 1 higher than given by
  // potentialCoefficients + 2 for the leading power, so its size + 2 (as it
  // gives coefficients for powers 2 to its size + 1)) which makes the
  // potential have a minimum when the auxiliary value is 1.0, which should
  // have been taken into account already in creating potentialCoefficients.
  inline void PathFieldsAndPotential::SetPotential(
                                 Eigen::VectorXd const& potentialCoefficients )
  {
    potentialApproximation = SimplePolynomial( potentialCoefficients,
                                               2,
                                               1 );
    std::vector< double >&
    potentialVector( potentialApproximation.CoefficientVector() );
    size_t const approximationDegree( potentialCoefficients.rows() + 1 );
    double finalCoefficientTimesMaxPower( 0.0 );
    for( size_t whichPower( 2 );
         whichPower < approximationDegree;
         ++whichPower )
    {
      finalCoefficientTimesMaxPower
      -= ( (double)whichPower * potentialVector[ whichPower ] );
    }
    potentialVector.back() = ( finalCoefficientTimesMaxPower
                               / (double)(approximationDegree) );
  }

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

} /* namespace VevaciousPlusPlus */
#endif /* PATHFIELDSANDPOTENTIAL_HPP_ */
