/*
 * BounceWithSplines.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{
  double const BounceWithSplines::maximumPowerOfNaturalExponent(
                           log( 0.5 * std::numeric_limits< double >::max() ) );
  double const
  BounceWithSplines::hBarInGigaElectronVoltSeconds( 6.58211928E-25 );
  double const BounceWithSplines::ageOfKnownUniverseInSeconds( 4.3E+17 );

  BounceWithSplines::BounceWithSplines( PotentialFunction& potentialFunction,
                                     TunnelingStrategy const tunnelingStrategy,
                                  double const survivalProbabilityThreshold ) :
    TunnelingCalculator( potentialFunction,
                         tunnelingStrategy,
                         survivalProbabilityThreshold ),
    potentialFunction( potentialFunction )
  {
    // This constructor is just an initialization list.
  }

  BounceWithSplines::~BounceWithSplines()
  {
    // This does nothing.
  }


  // This does something.
  void BounceWithSplines::CalculateTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // First we set all variables to their "not calculated" values.
    quantumSurvivalProbability = -1.0;
    quantumLifetimeInSeconds = -1.0;
    thermalSurvivalProbability = -1.0;
    dominantTemperatureInGigaElectronVolts = -1.0;

    if( tunnelingStrategy == NoTunneling )
    {
      std::cout
      << std::endl
      << "Not tunneling as tunneling strategy is \"NoTunneling\"";
      std::cout << std::endl;

      return;
    }
    PrepareCommonExtras();
    if( tunnelingStrategy == JustQuantum )
    {
      CalculateQuantumTunneling( falseVacuum,
                                 trueVacuum );
    }
    else if( tunnelingStrategy == JustThermal )
    {
      CalculateThermalTunneling( falseVacuum,
                                 trueVacuum );
    }
    else if( tunnelingStrategy == QuantumThenThermal )
    {
      CalculateQuantumTunneling( falseVacuum,
                                 trueVacuum );
      if( quantumSurvivalProbability > survivalProbabilityThreshold )
      {
        CalculateThermalTunneling( falseVacuum,
                                   trueVacuum );
      }
    }
    else if( tunnelingStrategy == ThermalThenQuantum )
    {
      CalculateThermalTunneling( falseVacuum,
                                 trueVacuum );
      if( thermalSurvivalProbability > survivalProbabilityThreshold )
      {
        CalculateQuantumTunneling( falseVacuum,
                                   trueVacuum );
      }
    }
    else
    {
      std::cout
      << std::endl
      << "No valid tunneling strategy was set, so treating it as"
      << " \"NoTunneling\"!";
      std::cout << std::endl;
    }
  }

} /* namespace VevaciousPlusPlus */
