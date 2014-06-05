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


  protected:
    SimplePolynomial potentialApproximation;
    std::vector< SimplePolynomial > fieldPath;
    std::vector< SimplePolynomial > pathTangent;
    double falseVacuumDepth;
    double trueVacuumDepth;
    bool nonZeroTemperature;
    double givenTemperature;
  };



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
