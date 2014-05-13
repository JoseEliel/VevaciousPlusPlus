/*
 * MinuitBounceActionMinimizer.hpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITBOUNCEACTIONMINIMIZER_HPP_
#define MINUITBOUNCEACTIONMINIMIZER_HPP_

#include "BounceWithSplines.hpp"
#include "../PotentialEvaluation.hpp"
#include "../PotentialMinimization.hpp"
#include "ModifiedBounceForMinuit.hpp"
#include "BubbleRadiusFromAuxiliary.hpp"

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


    // This should perform all relevant updates for the new SLHA data except
    // for propagating the push to the set of dependent SlhaUpdatePropagators.
    virtual void
    UpdateSelfForNewSlha( SlhaManager const& slhaManager );


  protected:
    // This should set quantumSurvivalProbability and quantumLifetimeInSeconds
    // appropriately.
    virtual void
    CalculateQuantumTunneling( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum );

    // This should set thermalSurvivalProbability and
    // dominantTemperatureInGigaElectronVolts appropriately.
    virtual void
    CalculateThermalTunneling( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum );
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINUITBOUNCEACTIONMINIMIZER_HPP_ */
