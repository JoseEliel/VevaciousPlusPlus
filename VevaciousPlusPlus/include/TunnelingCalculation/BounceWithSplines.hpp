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
    PotentialFunction const& potentialFunction;

    // This is a hook to allow for derived classes to prepare things common to
    // both quantum and thermal tunneling. By default, it does nothing.
    virtual void PrepareCommonExtras(){}

    // This should set quantumSurvivalProbability and quantumLifetimeInSeconds
    // appropriately.
    virtual void
    CalculateQuantumTunneling( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum ) = 0;

    // This should set thermalSurvivalProbability and
    // dominantTemperatureInGigaElectronVolts appropriately.
    virtual void
    CalculateThermalTunneling( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum ) = 0;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEWITHSPLINES_HPP_ */
