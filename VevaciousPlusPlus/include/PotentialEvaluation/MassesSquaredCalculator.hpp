/*
 * MassesSquaredCalculator.hpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MASSESSQUAREDCALCULATOR_HPP_
#define MASSESSQUAREDCALCULATOR_HPP_

#include "CommonIncludes.hpp"

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
    MassesSquaredCalculator( MassesSquaredCalculator const& copySource );
    MassesSquaredCalculator();
    virtual ~MassesSquaredCalculator();


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

    // This returns the eigenvalues of the matrix, with all functionoids
    // evaluated at the natural exponent of logarithmOfScale.
    virtual std::vector< double >
    MassesSquared( std::vector< double > const& fieldConfiguration,
                   double const logarithmOfScale ) const
    { // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Remember to remove this function!";
    std::cout << std::endl;/**/
    throw std::runtime_error( "no longer supported!" );
    return std::vector< double >();}

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
  };

} /* namespace VevaciousPlusPlus */
#endif /* MASSESSQUAREDCALCULATOR_HPP_ */
