/*
 * MinuitBounceActionMinimizer.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  MinuitBounceActionMinimizer::MinuitBounceActionMinimizer(
                                          PotentialFunction& potentialFunction,
                                     TunnelingStrategy const tunnelingStrategy,
                                  double const survivalProbabilityThreshold ) :
    BounceWithSplines( potentialFunction,
                       tunnelingStrategy,
                       survivalProbabilityThreshold )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::MinuitBounceActionMinimizer( ... )";
    std::cout << std::endl;/**/
  }

  MinuitBounceActionMinimizer::~MinuitBounceActionMinimizer()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::~MinuitBounceActionMinimizer()";
    std::cout << std::endl;/**/
  }


  // This doesn't do anything here.
  void MinuitBounceActionMinimizer::UpdateSelfForNewSlha(
                                               SlhaManager const& slhaManager )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::UpdateSelfForNewSlha(...)";
    std::cout << std::endl;/**/
  }

  // This should set quantumSurvivalProbability and quantumLifetimeInSeconds
  // appropriately.
  void MinuitBounceActionMinimizer::CalculateQuantumTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::CalculateQuantumTunneling(...)";
    std::cout << std::endl;/**/

    criticalRatherThanEvaporation = false;
    evaporationMinimum = falseVacuum;
    double evaporationTemperature( CriticalOrEvaporationTemperature(
                potentialFunction( potentialFunction.FieldValuesOrigin() ) ) );
    double
    appropriateScale( sqrt( potentialFunction.ScaleSquaredRelevantToTunneling(
                                                                   falseVacuum,
                                                              trueVacuum ) ) );
    BubbleRadiusFromAuxiliary bubbleRadiusFromAuxiliary( appropriateScale );
    ModifiedBounceForMinuit modifiedBounceForMinuit( potentialFunction,
                                             TunnelPathResolution( falseVacuum,
                                                                   trueVacuum,
                                                                   10 ),
                                                     falseVacuum,
                                                     evaporationTemperature,
                                                   bubbleRadiusFromAuxiliary );
    unsigned int const
    numberOfFields( potentialFunction.NumberOfFieldVariables() );
    std::vector< double > splineCoefficients( ( ( 2 * numberOfFields ) + 1 ),
                                              0.0 );
    for( unsigned int splineIndex( 0 );
         splineIndex < numberOfFields;
         ++splineIndex )
    {
      splineCoefficients[ splineIndex ]
      = ( trueVacuum.FieldConfiguration()[ splineIndex  ]
          - falseVacuum.FieldConfiguration()[ splineIndex ] );
    }
    double straightPathBounce( modifiedBounceForMinuit( splineCoefficients ) );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "straightPathBounce = " << straightPathBounce;
    std::cout << std::endl;/**/

    for( unsigned int splineIndex( 0 );
         splineIndex < numberOfFields;
         ++splineIndex )
    {
      splineCoefficients[ splineIndex ] = 0.0;
      splineCoefficients[ splineIndex + numberOfFields ]
      = ( trueVacuum.FieldConfiguration()[ splineIndex  ]
          - falseVacuum.FieldConfiguration()[ splineIndex ] );
    }
    straightPathBounce = modifiedBounceForMinuit( splineCoefficients );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "straightPathBounce = " << straightPathBounce;
    std::cout << std::endl;/**/
  }

  // This should set thermalSurvivalProbability and
  // dominantTemperatureInGigaElectronVolts appropriately.
  void MinuitBounceActionMinimizer::CalculateThermalTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::CalculateThermalTunneling(...)";
    std::cout << std::endl;/**/
  }
} /* namespace VevaciousPlusPlus */
