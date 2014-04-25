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


  // This should try to find the most accurate survival probability for
  // falseVacuum to have survived as long as the age of the known Universe
  // including the time at non-negligible temperatures, depending on
  // tunnelingStrategy. It should set quantumSurvivalProbability,
  // quantumLifetimeInSeconds, thermalSurvivalProbability, and
  // dominantTemperatureInGigaElectronVolts appropriately. Each of these
  // which is not calculated by the strategy should be left with negative
  // values.
  void MinuitBounceActionMinimizer::CalculateTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::CalculateTunneling()";
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
