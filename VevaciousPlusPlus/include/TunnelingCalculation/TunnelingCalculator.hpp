/*
 * TunnelingCalculator.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef TUNNELINGCALCULATOR_HPP_
#define TUNNELINGCALCULATOR_HPP_

#include "../CommonIncludes.hpp"
#include "../PotentialMinimization/PotentialMinimum.hpp"
#include "../PotentialEvaluation/SlhaUpdatePropagator.hpp"

namespace VevaciousPlusPlus
{
  class TunnelingCalculator : public SlhaUpdatePropagator
  {
  public:
    enum TunnelingStrategy
    {
      NoTunneling,
      JustQuantum,
      JustThermal,
      QuantumThenThermal,
      ThermalThenQuantum,
      NotSet
    };

    TunnelingCalculator( SlhaUpdatePropagator& slhaUpdatePropagator,
                         TunnelingStrategy const tunnelingStrategy,
                         double const survivalProbabilityThreshold );
    virtual
    ~TunnelingCalculator();


    // This should try to find the most accurate survival probability for
    // falseVacuum to have survived as long as the age of the known Universe
    // including the time at non-negligible temperatures, depending on
    // tunnelingStrategy. It should set quantumSurvivalProbability,
    // quantumLifetimeInSeconds, thermalSurvivalProbability, and
    // dominantTemperatureInGigaElectronVolts appropriately. Each of these
    // which is not calculated by the strategy should be left with negative
    // values.
    virtual void CalculateTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum ) = 0;


  protected:
    TunnelingStrategy const tunnelingStrategy;
    double quantumSurvivalProbability;
    double quantumLifetimeInSeconds;
    double thermalSurvivalProbability;
    double dominantTemperatureInGigaElectronVolts;
    double const survivalProbabilityThreshold;
  };

} /* namespace VevaciousPlusPlus */
#endif /* TUNNELINGCALCULATOR_HPP_ */
