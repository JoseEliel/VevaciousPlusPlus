/*
 * MinuitPotentialMinimizer.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/GradientBasedMinimization/MinuitPotentialMinimizer.hpp"

namespace VevaciousPlusPlus
{

  MinuitPotentialMinimizer::MinuitPotentialMinimizer(
                                    PotentialFunction const& potentialFunction,
                                                    double const errorFraction,
                                                    double const errorMinimum ,
                                          unsigned int const minuitStrategy ) :
    GradientMinimizer( potentialFunction ),
    minimizationFunction( potentialFunction ),
    errorFraction( errorFraction ),
    errorMinimum( errorMinimum ),
    minuitStrategy( minuitStrategy )
  {
    // This constructor is just an initialization list.
  }

  MinuitPotentialMinimizer::~MinuitPotentialMinimizer()
  {
    // This does nothing.
  }


  // This sets up a ROOT::Minuit2::MnMigrad instance and runs its operator().
  // The initial step sizes are set to be the values of startingPoint
  // multiplied by errorFraction, absolute values taken. Any step size less
  // than errorMinimum is set to errorMinimum.
  ROOT::Minuit2::FunctionMinimum MinuitPotentialMinimizer::RunMigrad(
                                    std::vector< double > const& startingPoint,
                                                  double givenTolerance ) const
  {
    std::vector< double > initialStepSizes( startingPoint.size() );
    double stepSize( 0.0 );
    for( size_t vectorIndex( 0 );
         vectorIndex < startingPoint.size();
         ++vectorIndex )
    {
      stepSize = ( errorFraction * startingPoint[ vectorIndex ] );
      if( stepSize < 0.0 )
      {
        stepSize = -stepSize;
      }
      if( stepSize < errorMinimum )
      {
        stepSize = errorMinimum;
      }
      initialStepSizes[ vectorIndex ] = stepSize;
    }
    if( givenTolerance <= 0.0 )
    {
      givenTolerance
      = std::max( errorMinimum,
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
