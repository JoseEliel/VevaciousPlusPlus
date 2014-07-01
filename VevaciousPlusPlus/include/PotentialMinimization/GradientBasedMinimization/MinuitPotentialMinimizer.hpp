/*
 * MinuitPotentialMinimizer.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITPOTENTIALMINIMIZER_HPP_
#define MINUITPOTENTIALMINIMIZER_HPP_

#include "CommonIncludes.hpp"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "../GradientMinimizer.hpp"
#include "MinuitMinimum.hpp"
#include "PotentialForMinuit.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"

namespace VevaciousPlusPlus
{

  class MinuitPotentialMinimizer : public GradientMinimizer
  {
  public:
    MinuitPotentialMinimizer( PotentialFunction const& potentialFunction,
                              std::string const& xmlArguments );
    MinuitPotentialMinimizer( PotentialFunction const& potentialFunction,
                              double const errorFraction = 0.1,
                              double const errorMinimum = 1.0,
                              unsigned int const minuitStrategy = 1 );
    virtual
    ~MinuitPotentialMinimizer();


    // This performs a Minuit2 migrad() minimization but puts the result in the
    // less cumbersome class PotentialMinimum instead of returning just a
    // ROOT::Minuit2::FunctionMinimum.
    virtual PotentialMinimum
    operator()( std::vector< double > const& startingPoint ) const
    { return PotentialMinimum( MinuitMinimum( startingPoint.size(),
                                              RunMigrad( startingPoint ) ) ); }

    // This sets up a ROOT::Minuit2::MnMigrad instance and runs its operator().
    // The initial step sizes are set to be the values of startingPoint
    // multiplied by errorFraction, absolute values taken. Any step size less
    // than errorMinimum is set to errorMinimum.
    ROOT::Minuit2::FunctionMinimum
    RunMigrad( std::vector< double > const& startingPoint,
               double givenTolerance = -1.0 ) const;


  protected:
    ROOT::Minuit2::FCNBase& minimizationFunction;
    double errorFraction;
    double errorMinimum;
    unsigned int minuitStrategy;
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINUITPOTENTIALMINIMIZER_HPP_ */
