/*
 * MinuitManager.hpp
 *
 *  Created on: Apr 9, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITMANAGER_HPP_
#define MINUITMANAGER_HPP_

#include "../../StandardIncludes.hpp"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "../PotentialMinimum.hpp"

namespace VevaciousPlusPlus
{
  // This is just a little wrapper class around ROOT::Minuit2::MnMigrad to ease
  // setting up the input.
  class MinuitManager
  {
  public:
    MinuitManager( ROOT::Minuit2::FCNBase& minimizationFunction,
                   double const errorFraction = 0.1,
                   double const errorMinimum = 1.0,
                   unsigned int const minuitStrategy = 1 );
    virtual
    ~MinuitManager();


    // This performs a Minuit2 MIGRAD minimization but puts the result in the
    // less cumbersome class MinuitMinimum instead of returning just a
    // ROOT::Minuit2::FunctionMinimum.
    MinuitMinimum
    operator()( std::vector< double > const& startingPoint,
                double givenTolerance = -1.0 ) const
    { return MinuitMinimum( startingPoint.size(),
                            RunMigrad( startingPoint,
                                       givenTolerance ) ); }

    // This sets up a ROOT::Minuit2::MnMigrad instance and runs its operator().
    // The initial step sizes are set to be the values of startingPoint
    // multiplied by errorFraction, absolute values taken. Any step size less
    // than errorMinimum is set to errorMinimum.
    ROOT::Minuit2::FunctionMinimum
    RunMigrad( std::vector< double > const& startingPoint,
               double givenTolerance = -1.0 ) const;


  protected:
    ROOT::Minuit2::FCNBase& minimizationFunction;
    double const errorFraction;
    double const errorMinimum;
    unsigned int const minuitStrategy;
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINUITMANAGER_HPP_ */
