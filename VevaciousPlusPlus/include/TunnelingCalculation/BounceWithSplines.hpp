/*
 * BounceWithSplines.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEWITHSPLINES_HPP_
#define BOUNCEWITHSPLINES_HPP_

#include "../CommonIncludes.hpp"
#include "../PotentialEvaluation.hpp"
#include "../PotentialMinimization.hpp"
#include "TunnelingCalculator.hpp"

namespace VevaciousPlusPlus
{
  class BounceWithSplines : public TunnelingCalculator
  {
  public:
    BounceWithSplines( PotentialFunction& potentialFunction,
                       TunnelingStrategy const tunnelingStrategy,
                       double const survivalProbabilityThreshold );
    virtual
    ~BounceWithSplines();


    // This decides what virtual tunneling calculation functions to call based
    // on tunnelingStrategy.
    virtual void CalculateTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum );


  protected:
    static double const maximumPowerOfNaturalExponent;
    static double const hBarInGigaElectronVoltSeconds;
    static double const ageOfKnownUniverseInSeconds;
    static double const ageOfKnownUniverseInInverseGigaElectronVolts;
    static double const fourVolumeOfKnownUniverseOverGevFourth;
    static double const lnOfThermalIntegrationFactor;
    PotentialFunction const& potentialFunction;
    PotentialForMinuit potentialForMinuit;
    MinuitManager thermalMinimumMinuit;
    PotentialMinimum evaporationMinimum;
    PotentialMinimum criticalMinimum;
    bool criticalRatherThanEvaporation;
    double thresholdSeparationSquared;

    // This is a hook to allow for derived classes to prepare things common to
    // both quantum and thermal tunneling. By default, it does nothing.
    virtual void PrepareCommonExtras(){}

    // This should return either the dimensionless bounce action integrated
    // over four dimensions (for zero temperature) or the dimensionful bounce
    // action integrated over three dimensions (for non-zero temperature) for
    // tunneling from falseVacuum to trueVacuum at temperature
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
    // vacuum.
    size_t TunnelPathResolution( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 size_t const resolutionOfDsbVacuum ) const
    { return ( resolutionOfDsbVacuum
              * (size_t)(sqrt( trueVacuum.SquareDistanceTo( falseVacuum )
                               / falseVacuum.LengthSquared() ) ) ); }
  };




  // This returns true if the temperature is below that at which tunneling
  // from the field origin to criticalMinimum becomes impossible.
  inline bool BounceWithSplines::BelowCriticalTemperature(
                                                double const temperatureGuess )
  {
    potentialForMinuit.SetTemperature( temperatureGuess );
    return ( potentialFunction( thermalMinimumMinuit(
                       criticalMinimum.FieldConfiguration() ).VariableValues(),
                                temperatureGuess )
             < potentialFunction( potentialFunction.FieldValuesOrigin(),
                                  temperatureGuess ) );
  }

  // This returns the result of BelowCriticalTemperature( temperatureGuess )
  // if criticalRatherThanEvaporation is true, otherwise the result of
  // BelowEvaporationTemperature( temperatureGuess ).
  inline bool BounceWithSplines::BelowCriticalOrEvaporation(
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

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEWITHSPLINES_HPP_ */
