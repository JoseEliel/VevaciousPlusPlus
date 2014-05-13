/*
 * ModifiedBounceForMinuit.hpp
 *
 *  Created on: May 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MODIFIEDBOUNCEFORMINUIT_HPP_
#define MODIFIEDBOUNCEFORMINUIT_HPP_

#include "../CommonIncludes.hpp"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "../PotentialEvaluation.hpp"
#include "BubbleRadiusFromAuxiliary.hpp"

namespace VevaciousPlusPlus
{

  class ModifiedBounceForMinuit : public ROOT::Minuit2::FCNBase
  {
  public:
    ModifiedBounceForMinuit( PotentialFunction const& potentialFunction,
                             unsigned int const numberOfPathNodes,
                             PotentialMinimum const& falseVacuum,
                             double const dsbEvaporationTemperature,
                  BubbleRadiusFromAuxiliary const& bubbleRadiusFromAuxiliary );
    virtual
    ~ModifiedBounceForMinuit();


    // This implements operator() for FCNBase, the function that MINUIT will
    // minimize. The values of splineCoefficients should be sets of n
    // coefficients for polynomials for each of the n fields, plus the final
    // element of splineCoefficients should be the temperature.
    virtual double
    operator()( std::vector< double > const& splineCoefficients ) const;

    // This implements Up() for FCNBase just to stick to a basic value.
    virtual double Up() const { return 1.0; }

    // This converts the spline coefficients into numberOfPathIntervals field
    // configurations and returns them as a vector of vectors.
    std::vector< std::vector< double > >
    PathFromSplines( std::vector< double > const& splineCoefficients ) const;


  protected:
    PotentialFunction const& potentialFunction;
    unsigned int numberOfFields;
    unsigned int numberOfPathIntervals;
    PotentialMinimum const& falseVacuum;
    double const falseVacuumEvaporationTemperature;
    BubbleRadiusFromAuxiliary const& bubbleRadiusFromAuxiliary;

    // This converts the spline coefficients into a field configuration and
    // puts it into fieldConfiguration.
    void ConfigurationFromSplines( std::vector< double >& fieldConfiguration,
                               std::vector< double > const& splineCoefficients,
                                   double const auxiliaryValue ) const;
  };




  // This converts the spline coefficients into differences from a starting
  // configuration given by fieldConfiguration and adds them to the values of
  // fieldConfiguration.
  inline void ModifiedBounceForMinuit::ConfigurationFromSplines(
                                     std::vector< double >& fieldConfiguration,
                               std::vector< double > const& splineCoefficients,
                                            double const auxiliaryValue ) const
  {
    // The last element of splineCoefficients is a temperature, so should not
    // be used.
    for( unsigned int splineIndex( 0 );
         splineIndex < ( splineCoefficients.size() - 1 );
         ++splineIndex )
    {
      // This relies heavily on int division.
      fieldConfiguration[ ( splineIndex % numberOfFields ) ]
      += ( splineCoefficients[ splineIndex ]
           * pow( auxiliaryValue,
                  ( ( splineIndex / numberOfFields ) + 1 ) ) );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* MODIFIEDBOUNCEFORMINUIT_HPP_ */
