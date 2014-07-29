/*
 * BounceActionTunneler.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEACTIONTUNNELINGCALCULATOR_HPP_
#define BOUNCEACTIONTUNNELINGCALCULATOR_HPP_

#include "CommonIncludes.hpp"
#include "TunnelingCalculator.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "SlhaManagement/SlhaManager.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitPotentialMinimizer.hpp"

namespace VevaciousPlusPlus
{

  class BounceActionTunneler : public TunnelingCalculator
  {
  public:
    BounceActionTunneler( PotentialFunction& potentialFunction,
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                          double const survivalProbabilityThreshold,
                          size_t const temperatureAccuracy,
                          size_t const evaporationResolution );
    virtual ~BounceActionTunneler();


    // This decides what virtual tunneling calculation functions to call based
    // on tunnelingStrategy.
    virtual void CalculateTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum );


    // This doesn't do anything here.
    virtual void UpdateSelfForNewSlha( SlhaManager const& slhaManager ){}


  protected:
    static double const maximumPowerOfNaturalExponent;
    static double const hBarInGigaElectronVoltSeconds;
    static double const ageOfKnownUniverseInSeconds;
    static double const ageOfKnownUniverseInInverseGigaElectronVolts;
    static double const fourVolumeOfKnownUniverseOverGevFourth;
    static double const lnOfThermalIntegrationFactor;

    PotentialFunction const& potentialFunction;
    size_t const temperatureAccuracy;
    size_t const evaporationResolution;
    MinuitPotentialMinimizer thermalPotentialMinimizer;
    PotentialMinimum evaporationMinimum;
    PotentialMinimum criticalMinimum;
    bool criticalRatherThanEvaporation;

    // This is a hook to allow for derived classes to prepare things common to
    // both quantum and thermal tunneling. By default, it does nothing.
    virtual void PrepareCommonExtras(){}

    // This should return either the dimensionless bounce action integrated
    // over four dimensions (for zero temperature) or the dimensionful bounce
    // action integrated over three dimensions (for non-zero temperature) for
    // tunneling from falseVacuum to trueVacuum at temperature
    // tunnelingTemperature. The vacua are assumed to already be the minima at
    // tunnelingTemperature.
    virtual double BounceAction( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 double const tunnelingTemperature ) const = 0;

    // This sets quantumSurvivalProbability and quantumLifetimeInSeconds
    // appropriately.
    virtual void
    CalculateQuantumTunneling( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum );

    // This should set thermalSurvivalProbability and
    // dominantTemperatureInGigaElectronVolts appropriately.
    virtual void
    CalculateThermalTunneling( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum ) = 0;

    // This calculates the temperature at which either tunneling from
    // givenVacuum to the field origin becomes impossible if
    // criticalRatherThanEvaporation is true or the temperature at which
    // givenVacuum evaporates if false.
    double
    CriticalOrEvaporationTemperature( double const potentialAtOrigin );

    // This returns true if the temperature is below that at which tunneling
    // from the field origin to criticalMinimum becomes impossible.
    bool BelowCriticalTemperature( double const temperatureGuess );

    // This returns true if the temperature is below that at which
    // evaporationMinimum is no longer separated from the field origin by an
    // energy barrier.
    bool BelowEvaporationTemperature( double const temperatureGuess );

    // This returns the result of BelowCriticalTemperature( temperatureGuess )
    // if criticalRatherThanEvaporation is true, otherwise the result of
    // BelowEvaporationTemperature( temperatureGuess ).
    bool BelowCriticalOrEvaporation( double const temperatureGuess );

    // This returns a number of points which should be appropriate for
    // resolving the potential to the extent that there are
    // resolutionOfDsbVacuum points between the field origin and the false
    // vacuum (if the false vacuum is not the field origin: if it is, just
    // resolutionOfDsbVacuum is returned).
    size_t TunnelPathResolution( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 size_t const resolutionOfDsbVacuum ) const;

    // This ensures that thermalSurvivalProbability is set correctly from
    // logOfMinusLogOfThermalProbability.
    void SetThermalSurvivalProbability();
  };




  // This returns true if the temperature is below that at which tunneling
  // from the field origin to criticalMinimum becomes impossible.
  inline bool BounceActionTunneler::BelowCriticalTemperature(
                                                double const temperatureGuess )
  {
    thermalPotentialMinimizer.SetTemperature( temperatureGuess );
    return ( potentialFunction( thermalPotentialMinimizer(
                       criticalMinimum.FieldConfiguration() ).VariableValues(),
                                temperatureGuess )
             < potentialFunction( potentialFunction.FieldValuesOrigin(),
                                  temperatureGuess ) );
  }

  // This returns the result of BelowCriticalTemperature( temperatureGuess )
  // if criticalRatherThanEvaporation is true, otherwise the result of
  // BelowEvaporationTemperature( temperatureGuess ).
  inline bool BounceActionTunneler::BelowCriticalOrEvaporation(
                                                double const temperatureGuess )
  {
    if( criticalRatherThanEvaporation )
    {
      return BelowCriticalTemperature( temperatureGuess );
    }
    else
    {
      return BelowEvaporationTemperature( temperatureGuess );
    }
  }

  // This returns a number of points which should be appropriate for
  // resolving the potential to the extent that there are
  // resolutionOfDsbVacuum points between the field origin and the false
  // vacuum (if the false vacuum is not the field origin: if it is, just
  // resolutionOfDsbVacuum is returned).
  inline size_t BounceActionTunneler::TunnelPathResolution(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                     size_t const resolutionOfDsbVacuum ) const
  {
    double const falseVacuumLengthSquared( falseVacuum.LengthSquared() );
    if( falseVacuumLengthSquared > 0.0 )
    {
      return ( resolutionOfDsbVacuum
               * (size_t)(sqrt( trueVacuum.SquareDistanceTo( falseVacuum )
                                / falseVacuumLengthSquared ) ) );
    }
    else
    {
      return resolutionOfDsbVacuum;
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEACTIONTUNNELINGCALCULATOR_HPP_ */
