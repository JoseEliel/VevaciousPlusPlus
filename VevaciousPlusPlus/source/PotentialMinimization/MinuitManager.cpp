/*
 * MinuitManager.cpp
 *
 *  Created on: Apr 9, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  MinuitManager::MinuitManager( ROOT::Minuit2::FCNBase& minimizationFunction,
                                double const errorFraction,
                                double const errorMinimum,
                                unsigned int const minuitStrategy ) :
    minimizationFunction( minimizationFunction ),
    errorFraction( errorFraction ),
    errorMinimum( errorMinimum ),
    minuitStrategy( minuitStrategy )
  {
    // This constructor is just an initialization list.
  }

  MinuitManager::~MinuitManager()
  {
    // This does nothing.
  }


  // This sets up a ROOT::Minuit2::MnMigrad instance and runs its operator().
  // The initial step sizes are set to be the values of startingPoint
  // multiplied by errorFraction, absolute values taken. Any step size less
  // than errorMinimum is set to errorMinimum.
  ROOT::Minuit2::FunctionMinimum
  MinuitManager::RunMigrad( std::vector< double > const& startingPoint,
                             double givenTolerance ) const
  {
    std::vector< double > initialStepSizes( startingPoint.size() );
    double stepSize( 0.0 );
    for( unsigned int vectorIndex( 0 );
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
      = ( errorFraction * minimizationFunction( startingPoint ) );
    }
    if( givenTolerance <= 0.0 )
    {
      givenTolerance = errorMinimum;
    }
    ROOT::Minuit2::MnMigrad mnMigrad( minimizationFunction,
                                      startingPoint,
                                      initialStepSizes,
                                      minuitStrategy );
    return mnMigrad( 0,
                     givenTolerance );
  }

} /* namespace VevaciousPlusPlus */
