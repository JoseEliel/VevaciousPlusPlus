/*
 * TunnelingCalculator.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  TunnelingCalculator::TunnelingCalculator(
                                          PotentialFunction& potentialFunction,
                                     TunnelingStrategy const tunnelingStrategy,
                                  double const survivalProbabilityThreshold ) :
    potentialFunction( potentialFunction ),
    tunnelingStrategy( tunnelingStrategy ),
    quantumSurvivalProbability( -1.0 ),
    quantumLifetimeInSeconds( -1.0 ),
    thermalSurvivalProbability( -1.0 ),
    dominantTemperatureInGigaElectronVolts( -1.0 ),
    survivalProbabilityThreshold( survivalProbabilityThreshold )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "TunnelingCalculator::TunnelingCalculator()";
    std::cout << std::endl;/**/
  }

  TunnelingCalculator::~TunnelingCalculator()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
