/*
 * TunnelingCalculator.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/TunnelingCalculator.hpp"

namespace VevaciousPlusPlus
{

  TunnelingCalculator::TunnelingCalculator(
                                    SlhaUpdatePropagator& slhaUpdatePropagator,
                                     TunnelingStrategy const tunnelingStrategy,
                                  double const survivalProbabilityThreshold ) :
    SlhaUpdatePropagator( slhaUpdatePropagator ),
    tunnelingStrategy( tunnelingStrategy ),
    quantumSurvivalProbability( -1.0 ),
    logOfMinusLogOfQuantumProbability( -1.0E+100 ),
    quantumLifetimeInSeconds( -1.0 ),
    thermalSurvivalProbability( -1.0 ),
    logOfMinusLogOfThermalProbability( -1.0E+100 ),
    dominantTemperatureInGigaElectronVolts( -1.0 ),
    survivalProbabilityThreshold( survivalProbabilityThreshold )
  {
    // This constructor is just an initialization list.
  }

  TunnelingCalculator::~TunnelingCalculator()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
