/*
 * BounceActionTunneler.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEACTIONTUNNELINGCALCULATOR_HPP_
#define BOUNCEACTIONTUNNELINGCALCULATOR_HPP_

#include "TunnelingCalculator.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include <utility>
#include <iostream>
#include "PotentialMinimization/GradientBasedMinimization/MinuitPotentialMinimizer.hpp"
#include <cmath>
#include <limits>
#include "Utilities/WarningLogger.hpp"
#include <vector>

namespace VevaciousPlusPlus
{

  class BounceActionTunneler : public TunnelingCalculator
  {
  public:
    BounceActionTunneler(
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                          double const survivalProbabilityThreshold,
                          unsigned int const temperatureAccuracy,
                          double const vacuumSeparationFraction );
    virtual ~BounceActionTunneler();


    // This decides what virtual tunneling calculation functions to call based
    // on tunnelingStrategy.
    virtual void
    CalculateTunneling( PotentialFunction const& potentialFunction,
                        PotentialMinimum const& falseVacuum,
                        PotentialMinimum const& trueVacuum );


  protected:
    static double const maximumPowerOfNaturalExponent;
    // The maximum allowed temperature is the reduced Planck mass.
    static double const maximumAllowedTemperature;
    static double const hBarInGigaElectronVoltSeconds;
    static double const ageOfKnownUniverseInSeconds;
    static double const ageOfKnownUniverseInInverseGigaElectronVolts;
    static double const fourVolumeOfKnownUniverseOverGevFourth;
    static double const lnOfThermalIntegrationFactor;

    unsigned int const temperatureAccuracy;
    std::pair< double, double > rangeOfMaxTemperatureForOriginToFalse;
    std::pair< double, double > rangeOfMaxTemperatureForOriginToTrue;

    double const vacuumSeparationFractionSquared;


    // This is a hook to allow for derived classes to prepare things common to
    // both quantum and thermal tunneling. By default, it does nothing.
    virtual void
    PrepareCommonExtras( PotentialFunction const& potentialFunction )
    { ; /* By default, does nothing. */ }

    // This should return either the dimensionless bounce action integrated
    // over four dimensions (for zero temperature) or the dimensionful bounce
    // action integrated over three dimensions (for non-zero temperature) for
    // tunneling from falseVacuum to trueVacuum at temperature
    // tunnelingTemperature. The vacua are assumed to already be the minima at
    // tunnelingTemperature.
    virtual double BounceAction( PotentialFunction const& potentialFunction,
                                 PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 double const tunnelingTemperature ) = 0;

    // This sets quantumSurvivalProbability, quantumLifetimeInSeconds, and
    // logOfMinusLogOfQuantumProbability appropriately.
    virtual void
    CalculateQuantumTunneling( PotentialFunction const& potentialFunction,
                               PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum );

    // This sets thermalSurvivalProbability,
    // dominantTemperatureInGigaElectronVolts, and
    // logOfMinusLogOfThermalProbability appropriately.
    virtual void
    CalculateThermalTunneling( PotentialFunction const& potentialFunction,
                               PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum );

    // This sets rangeOfMaxTemperatureForOriginToFalse and
    // rangeOfMaxTemperatureForOriginToTrue to be pairs of temperatures which
    // are just above and just below the maximum temperatures for tunneling to
    // be possible from the origin to the false vacuum and true vacuum
    // respectively. The temperatures are capped at the Planck temperature.
    void
    SetUpMaximumTemperatureRanges( PotentialFunction const& potentialFunction,
                                   PotentialMinimum const& falseVacuum,
                                   PotentialMinimum const& trueVacuum,
                             double const potentialAtOriginAtZeroTemperature );

    // This sets rangeOfMaxTemperature to be the temperatures which are just
    // above and just below the maximum temperature for tunneling to be
    // possible from the origin to zeroTemperatureVacuum. The temperatures are
    // capped at the Planck temperature.
    void SetMaximumTunnelingTemperatureRange(
                                    PotentialFunction const& potentialFunction,
                            std::pair< double, double >& rangeOfMaxTemperature,
                                 PotentialMinimum const& zeroTemperatureVacuum,
                             double const potentialAtOriginAtZeroTemperature );

    // This should set thermalSurvivalProbability,
    // dominantTemperatureInGigaElectronVolts, and
    // logOfMinusLogOfThermalProbability appropriately, in the context of being
    // called after guarding against the case where the field origin is deeper
    // than the DSB vacuum at zero temperature, which would make cooling into
    // the DSB vacuum implausible.
    virtual void
    ContinueThermalTunneling( PotentialFunction const& potentialFunction,
                              PotentialMinimum const& falseVacuum,
                              PotentialMinimum const& trueVacuum,
                         double const potentialAtOriginAtZeroTemperature ) = 0;

    // This returns true if the temperature is below that at which tunneling
    // from the field origin to zeroTemperatureVacuum becomes impossible.
    bool BelowCriticalTemperature( PotentialFunction const& potentialFunction,
                                   double const temperatureGuess,
                               PotentialMinimum const& zeroTemperatureVacuum );

    // This returns a number of points which should be appropriate for
    // resolving the potential to the extent that there are
    // resolutionOfDsbVacuum points between the field origin and the false
    // vacuum (if the false vacuum is not the field origin: if it is, just
    // resolutionOfDsbVacuum is returned).
    unsigned int TunnelPathResolution( PotentialMinimum const& falseVacuum,
                                       PotentialMinimum const& trueVacuum,
                              unsigned int const resolutionOfDsbVacuum ) const;

    // This ensures that thermalSurvivalProbability is set correctly from
    // logOfMinusLogOfThermalProbability.
    void SetThermalSurvivalProbability();

    // This should return the A factor of dimension energy^4.
    virtual double SolitonicFactor( PotentialFunction const& potentialFunction,
                                    PotentialMinimum const& falseVacuum,
                                    PotentialMinimum const& trueVacuum );
  };




  // This sets rangeOfMaxTemperatureForOriginToFalse and
  // rangeOfMaxTemperatureForOriginToTrue to be pairs of temperatures which
  // are just above and just below the maximum temperatures for tunneling to
  // be possible from the origin to the false vacuum and true vacuum
  // respectively. The temperatures are capped at the Planck temperature.
  inline void BounceActionTunneler::SetUpMaximumTemperatureRanges(
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                              double const potentialAtOriginAtZeroTemperature )
  {
    std::cout << std::endl
    << "Looking for temperature at which tunneling from the field origin to"
    << " the false vacuum at "
    << falseVacuum.AsMathematica( potentialFunction.FieldNames() )
    << " becomes impossible." << std::endl;
    SetMaximumTunnelingTemperatureRange( potentialFunction,
                                         rangeOfMaxTemperatureForOriginToFalse,
                                         falseVacuum,
                                         potentialAtOriginAtZeroTemperature );
    std::cout << std::endl
    << "Looking for temperature at which tunneling from the field origin to"
    << " the true vacuum at "
    << trueVacuum.AsMathematica( potentialFunction.FieldNames() )
    << " becomes impossible." << std::endl;
    SetMaximumTunnelingTemperatureRange( potentialFunction,
                                         rangeOfMaxTemperatureForOriginToTrue,
                                         trueVacuum,
                                         potentialAtOriginAtZeroTemperature );
  }

  // This returns true if the temperature is below that at which tunneling
  // from the field origin to zeroTemperatureVacuum becomes impossible.
  inline bool BounceActionTunneler::BelowCriticalTemperature(
                                    PotentialFunction const& potentialFunction,
                                                 double const temperatureGuess,
                                PotentialMinimum const& zeroTemperatureVacuum )
  {
    double const thermalPotentialAtOrigin( potentialFunction(
                                         potentialFunction.FieldValuesOrigin(),
                                           temperatureGuess ) );
    if( potentialFunction( zeroTemperatureVacuum.FieldConfiguration(),
                           temperatureGuess  ) < thermalPotentialAtOrigin )
    {
      return true;
    }
    else
    {
      MinuitPotentialMinimizer thermalPotentialMinimizer( potentialFunction );
      thermalPotentialMinimizer.SetTemperature( temperatureGuess );
      return ( ( thermalPotentialMinimizer(
                   zeroTemperatureVacuum.FieldConfiguration() ).FunctionValue()
                 + thermalPotentialMinimizer.FunctionOffset() )
               < thermalPotentialAtOrigin );
    }
  }

  // This returns a number of points which should be appropriate for
  // resolving the potential to the extent that there are
  // resolutionOfDsbVacuum points between the field origin and the false
  // vacuum (if the false vacuum is not the field origin: if it is, just
  // resolutionOfDsbVacuum is returned).
  inline unsigned int BounceActionTunneler::TunnelPathResolution(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                               unsigned int const resolutionOfDsbVacuum ) const
  {
    double const falseVacuumLengthSquared( falseVacuum.LengthSquared() );
    if( falseVacuumLengthSquared > 0.0 )
    {
      return ( resolutionOfDsbVacuum
               * static_cast< unsigned int >(
                               sqrt( trueVacuum.SquareDistanceTo( falseVacuum )
                                     / falseVacuumLengthSquared ) ) );
    }
    else
    {
      return resolutionOfDsbVacuum;
    }
  }

  // This should return the A factor of dimension energy^4. This particular
  // implementation does the naive thing of just taking an appropriate energy
  // scale and returning it to the fourth power.
  inline double BounceActionTunneler::SolitonicFactor(
                                   PotentialFunction const& potentialFunction,
                                  PotentialMinimum const& falseVacuum,
                                  PotentialMinimum const& trueVacuum )
  {
    double const squareRootOfSolitonicFactor(
              potentialFunction.ScaleSquaredRelevantToTunneling( falseVacuum,
                                                                trueVacuum ) );
    return ( squareRootOfSolitonicFactor * squareRootOfSolitonicFactor );
  }

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEACTIONTUNNELINGCALCULATOR_HPP_ */
