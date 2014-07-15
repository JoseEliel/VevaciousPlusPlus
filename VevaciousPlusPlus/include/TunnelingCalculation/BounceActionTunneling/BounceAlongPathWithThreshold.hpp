/*
 * BounceAlongPathWithThreshold.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEALONGPATHWITHTHRESHOLD_HPP_
#define BOUNCEALONGPATHWITHTHRESHOLD_HPP_

#include "CommonIncludes.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "../BounceActionTunneler.hpp"
#include "BounceActionEvaluation/BouncePathFinder.hpp"
#include "BounceActionEvaluation/BounceActionCalculator.hpp"

namespace VevaciousPlusPlus
{
  // This class takes memory management ownership of the components given to
  // the constructor as pointers! It'd be nice to use std::unique_ptrs, but we
  // are stubbornly sticking to allowing non-C++11-compliant compilers.
  class BounceAlongPathWithThreshold : public BounceActionTunneler
  {
  public:
    BounceAlongPathWithThreshold( PotentialFunction& potentialFunction,
                                  BouncePathFinder* pathFinder,
                                  BounceActionCalculator* actionCalculator,
                                  std::string const& xmlArguments );
    virtual ~BounceAlongPathWithThreshold();


  protected:
    BouncePathFinder* pathFinder;
    BounceActionCalculator* actionCalculator;
    double actionThreshold;
    size_t thermalIntegrationResolution;


    // This returns either the dimensionless bounce action integrated over four
    // dimensions (for zero temperature) or the dimensionful bounce action
    // integrated over three dimensions (for non-zero temperature) for
    // tunneling from falseVacuum to trueVacuum at temperature
    // tunnelingTemperature, or an upper bound if the upper bound drops below
    // actionThreshold during the course of the calculation. The vacua are
    // assumed to already be the minima at tunnelingTemperature.
    virtual double BounceAction( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 double const tunnelingTemperature ) const;

    // This sets thermalSurvivalProbability by numerically integrating from the
    // critical temperature for tunneling to be possible down to T = 0 unless
    // the integral already passes a threshold, and sets
    // dominantTemperatureInGigaElectronVolts to be the temperature with the
    // lowest survival probability.
    virtual void
    CalculateThermalTunneling( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum );
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEALONGPATHWITHTHRESHOLD_HPP_ */
