/*
 * MinuitPotentialMinimizer.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITPOTENTIALMINIMIZER_HPP_
#define MINUITPOTENTIALMINIMIZER_HPP_

#include "PotentialMinimization/GradientMinimizer.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include <vector>
#include "MinuitWrappersAndHelpers/MinuitMinimum.hpp"
#include "Minuit2/FunctionMinimum.h"
#include "PotentialForMinuit.hpp"
#include "Minuit2/MnMigrad.h"
#include <cstddef>
#include <algorithm>
#include <cmath>


namespace VevaciousPlusPlus
{

  class MinuitPotentialMinimizer : public GradientMinimizer
  {
  public:
    MinuitPotentialMinimizer( PotentialFunction const& potentialFunction,
                              double const errorFraction = 0.1,
                              double const errorMinimum = 1.0,
                              unsigned int const minuitStrategy = 1 ) :
      GradientMinimizer( potentialFunction ),
      minimizationFunction( potentialFunction ),
      errorFraction( errorFraction ),
      errorMinimum( errorMinimum ),
      minuitStrategy( minuitStrategy ) {}

    virtual ~MinuitPotentialMinimizer() {}


    // This performs a Minuit2 migrad() minimization but puts the result in the
    // less cumbersome class PotentialMinimum instead of returning just a
    // ROOT::Minuit2::FunctionMinimum.
    virtual PotentialMinimum
    operator()( std::vector< double > const& startingPoint ) const
    { return PotentialMinimum( MinuitMinimum( startingPoint.size(),
                                              RunMigrad( startingPoint ) ) ); }

    // This ensures that the minimizations are calculated at the given
    // temperature.
    virtual void SetTemperature( double const minimizationTemperature )
    { minimizationFunction.SetTemperature( minimizationTemperature ); }

    // This sets up a ROOT::Minuit2::MnMigrad instance and runs its operator().
    // The initial step sizes are set to be the values of startingPoint
    // multiplied by errorFraction, absolute values taken. Any step size less
    // than errorMinimum is set to errorMinimum.
    ROOT::Minuit2::FunctionMinimum
    RunMigrad( std::vector< double > const& startingPoint,
               double givenTolerance = -1.0 ) const;

    // This returns the value of the potential at the field origin and at the
    // current temperature, which is subtracted from the potential by
    // minimizationFunction when evaluating the potential for Minuit2.
    double FunctionOffset() const
    { return minimizationFunction.FunctionAtOrigin(); }


  protected:
    PotentialForMinuit minimizationFunction;
    double const errorFraction;
    double const errorMinimum;
    unsigned int const minuitStrategy;
  };





  // This sets up a ROOT::Minuit2::MnMigrad instance and runs its operator().
  // The initial step sizes are set to be the values of startingPoint
  // multiplied by errorFraction, absolute values taken. Any step size less
  // than errorMinimum is set to errorMinimum.
  inline ROOT::Minuit2::FunctionMinimum MinuitPotentialMinimizer::RunMigrad(
                                    std::vector< double > const& startingPoint,
                                                  double givenTolerance ) const
  {
    std::vector< double > initialStepSizes( startingPoint.size(),
                                            errorMinimum );
    for( size_t vectorIndex( 0 );
         vectorIndex < startingPoint.size();
         ++vectorIndex )
    {
      initialStepSizes[ vectorIndex ] = std::max( errorMinimum,
                        fabs( errorFraction * startingPoint[ vectorIndex ] ) );
    }
    if( givenTolerance <= 0.0 )
    {
      givenTolerance = std::max( errorMinimum,
                  ( errorFraction * minimizationFunction( startingPoint ) ) );
    }
    ROOT::Minuit2::MnMigrad mnMigrad( minimizationFunction,
                                      startingPoint,
                                      initialStepSizes,
                                      minuitStrategy );
    return mnMigrad( 0,
                     givenTolerance );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITPOTENTIALMINIMIZER_HPP_ */
