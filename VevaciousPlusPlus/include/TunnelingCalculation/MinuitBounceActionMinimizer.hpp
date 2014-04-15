/*
 * MinuitBounceActionMinimizer.hpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITBOUNCEACTIONMINIMIZER_HPP_
#define MINUITBOUNCEACTIONMINIMIZER_HPP_

namespace VevaciousPlusPlus
{

  class MinuitBounceActionMinimizer : public BounceWithSplines
  {
  public:
    MinuitBounceActionMinimizer( PotentialFunction& potentialFunction,
                                 TunnelingStrategy const tunnelingStrategy,
                                 double const survivalProbabilityThreshold );
    virtual
    ~MinuitBounceActionMinimizer();


    // This should try to find the most accurate survival probability for
    // falseVacuum to have survived as long as the age of the known Universe
    // including the time at non-negligible temperatures, depending on
    // tunnelingStrategy. It should set quantumSurvivalProbability,
    // quantumLifetimeInSeconds, thermalSurvivalProbability, and
    // dominantTemperatureInGigaElectronVolts appropriately. Each of these
    // which is not calculated by the strategy should be left with negative
    // values.
    virtual void CalculateTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum );
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINUITBOUNCEACTIONMINIMIZER_HPP_ */
