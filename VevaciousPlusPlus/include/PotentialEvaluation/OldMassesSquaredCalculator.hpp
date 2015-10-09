/*
 * MassesSquaredCalculator.hpp
 *
 *  Created on: Apr 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef OLDMASSESSQUAREDCALCULATOR_HPP_
#define OLDMASSESSQUAREDCALCULATOR_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class OldMassesSquaredCalculator
  {
  public:
    enum SpinType
    {
      scalarBoson,
      weylFermion,
      gaugeBoson,
      notSet
    };

    OldMassesSquaredCalculator(
                    std::map< std::string, std::string > const& attributeMap );
    OldMassesSquaredCalculator( OldMassesSquaredCalculator const& copySource );
    OldMassesSquaredCalculator();
    virtual ~OldMassesSquaredCalculator();


    // This should return the masses-squared, with all functionoids evaluated
    // at the last scale which was used to update them.
    virtual std::vector< double >
    MassesSquared( std::vector< double > const& fieldConfiguration ) const = 0;

    // This should return the masses-squared, with all functionoids evaluated
    // at the natural exponent of logarithmOfScale.
    virtual std::vector< double >
    MassesSquared( std::vector< double > const& fieldConfiguration,
                   double const logarithmOfScale ) const = 0;

    // This returns the number of identical copies of this mass-squared matrix
    // that the model has.
    double MultiplicityFactor() const{ return multiplicityFactor; }

    SpinType GetSpinType() const{ return spinType; }


  protected:
    double multiplicityFactor;
    SpinType spinType;
  };

} /* namespace VevaciousPlusPlus */
#endif /* OLDMASSESSQUAREDCALCULATOR_HPP_ */
