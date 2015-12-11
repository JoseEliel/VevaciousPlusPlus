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
