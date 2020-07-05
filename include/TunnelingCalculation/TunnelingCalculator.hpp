/*
 * TunnelingCalculator.hpp
 *
 *  Created on: Feb 25, 2014
 *      Authors: Ben O'Leary (benjamin.oleary@gmail.com)
 *               José Eliel Camargo-Molina (elielcamargomolina@gmail.com)
 */

#ifndef TUNNELINGCALCULATOR_HPP_
#define TUNNELINGCALCULATOR_HPP_

#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"

namespace VevaciousPlusPlus
{

  class TunnelingCalculator
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

    TunnelingCalculator( TunnelingStrategy const tunnelingStrategy = NotSet,
                         double const survivalProbabilityThreshold = 0.01 ) :
      tunnelingStrategy( tunnelingStrategy ),
      quantumSurvivalProbability( -1.0 ),
      logOfMinusLogOfQuantumProbability( -1.0E+100 ),
      quantumLifetimeInSeconds( -1.0 ),
      thermalSurvivalProbability( -1.0 ),
      logOfMinusLogOfThermalProbability( -1.0E+100 ),
      partialThermalDecayWidth(-1),
      dominantTemperatureInGigaElectronVolts( -1.0 ),
      survivalProbabilityThreshold( survivalProbabilityThreshold ),
      thresholdAndActions(0),
      thermalThresholdAndActions(0) {}

    virtual ~TunnelingCalculator() {}


    // This should try to find the most accurate survival probability for
    // falseVacuum to have survived as long as the age of the known Universe
    // including the time at non-negligible temperatures, depending on
    // tunnelingStrategy. It should set quantumSurvivalProbability,
    // logOfMinusLogOfQuantumProbability, quantumLifetimeInSeconds,
    // thermalSurvivalProbability, logOfMinusLogOfThermalProbability, and
    // dominantTemperatureInGigaElectronVolts appropriately. Each of these
    // which is not calculated by the strategy should be left with negative
    // values.
    virtual void
    CalculateTunneling( PotentialFunction const& potentialFunction,
                        PotentialMinimum const& falseVacuum,
                        PotentialMinimum const& trueVacuum ) = 0;

    double QuantumSurvivalProbability() const
    { return quantumSurvivalProbability; }

    double LogOfMinusLogOfQuantumProbability() const
    { return logOfMinusLogOfQuantumProbability; }

    double QuantumLifetimeInSeconds() const
    { return quantumLifetimeInSeconds; }

    double ThermalSurvivalProbability() const
    { return thermalSurvivalProbability; }

    double LogOfMinusLogOfThermalProbability() const
    { return logOfMinusLogOfThermalProbability; }

    double DominantTemperatureInGigaElectronVolts() const
    { return dominantTemperatureInGigaElectronVolts; }

    double SurvivalProbabilityThreshold() const
    { return survivalProbabilityThreshold; }

    double PartialThermalDecayWidth() const
    { return partialThermalDecayWidth; }

    // This returns a vector with the bounce action for the straight path
    // and the best action found by each used pathfinder.

    std::vector< double > GetThresholdAndActions() const
    { return thresholdAndActions; }

    std::vector< double > GetThermalThresholdAndActions() const
    { return thermalThresholdAndActions;}


  protected:
    TunnelingStrategy tunnelingStrategy;
    double quantumSurvivalProbability;
    double logOfMinusLogOfQuantumProbability;
    double quantumLifetimeInSeconds;
    double thermalSurvivalProbability;
    double logOfMinusLogOfThermalProbability;
    double dominantTemperatureInGigaElectronVolts;
    double survivalProbabilityThreshold;
    double partialThermalDecayWidth;
    // The vectors below will hold the action threshold in their 0 component
    // the straight path action in their 1 component
    // and subsequently the best action for each path finder
    std::vector< double > thresholdAndActions;
    std::vector< double > thermalThresholdAndActions;
  };

} /* namespace VevaciousPlusPlus */
#endif /* TUNNELINGCALCULATOR_HPP_ */
