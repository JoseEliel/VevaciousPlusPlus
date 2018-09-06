/*
 * BounceAlongPathWithThreshold.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEALONGPATHWITHTHRESHOLD_HPP_
#define BOUNCEALONGPATHWITHTHRESHOLD_HPP_

#include "TunnelingCalculation/BounceActionTunneler.hpp"
#include <vector>
#include "BounceActionEvaluation/BouncePathFinder.hpp"
#include "BounceActionEvaluation/BounceActionCalculator.hpp"
#include "TunnelingCalculation/TunnelingCalculator.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include <cmath>
#include <cstddef>
#include "PotentialMinimization/GradientBasedMinimization/MinuitPotentialMinimizer.hpp"
#include <iostream>
#include "Utilities/WarningLogger.hpp"
#include "BounceActionEvaluation/PathParameterization/TunnelPath.hpp"
#include "BounceActionEvaluation/PathParameterization/LinearSplineThroughNodes.hpp"
#include "BounceActionEvaluation/SplinePotential.hpp"
#include "BounceActionEvaluation/BubbleProfile.hpp"

namespace VevaciousPlusPlus
{
  class BounceAlongPathWithThreshold : public BounceActionTunneler
  {
  public:
    BounceAlongPathWithThreshold(
                           std::vector< std::shared_ptr<BouncePathFinder> > const& pathFinders,
                                std::unique_ptr<BounceActionCalculator> actionCalculator,
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                                  double const survivalProbabilityThreshold,
                               unsigned int const thermalIntegrationResolution,
                                  unsigned int const temperatureAccuracy,
                                  unsigned int const pathPotentialResolution,
                                  double const vacuumSeparationFraction );
    virtual ~BounceAlongPathWithThreshold();


  protected:
    std::vector< std::shared_ptr<BouncePathFinder> > pathFinders;
    std::unique_ptr<BounceActionCalculator> actionCalculator;
    unsigned int thermalIntegrationResolution;
    unsigned int const pathPotentialResolution;

    // This returns either the dimensionless bounce action integrated over four
    // dimensions (for zero temperature) or the dimensionful bounce action
    // integrated over three dimensions (for non-zero temperature) for
    // tunneling from falseVacuum to trueVacuum at temperature
    // tunnelingTemperature, or an upper bound if the upper bound drops below
    // actionThreshold during the course of the calculation. The vacua are
    // assumed to already be the minima at tunnelingTemperature.
    virtual double BounceAction( PotentialFunction const& potentialFunction,
                                 PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 double const tunnelingTemperature );

    // This sets thermalSurvivalProbability by numerically integrating up to
    // the critical temperature for tunneling to be possible from T = 0 unless
    // the integral already passes a threshold, and sets
    // dominantTemperatureInGigaElectronVolts to be the temperature with the
    // lowest survival probability.
    virtual void
    ContinueThermalTunneling( PotentialFunction const& potentialFunction,
                              PotentialMinimum const& falseVacuum,
                              PotentialMinimum const& trueVacuum,
                             double const potentialAtOriginAtZeroTemperature );

    // This returns either the dimensionless bounce action integrated over four
    // dimensions (for zero temperature) or the dimensionful bounce action
    // integrated over three dimensions (for non-zero temperature) for
    // tunneling from falseVacuum to trueVacuum at temperature
    // tunnelingTemperature, or an upper bound if the upper bound drops below
    // actionThreshold during the course of the calculation. The vacua are
    // assumed to already be the minima at tunnelingTemperature.
    double BoundedBounceAction( PotentialFunction const& potentialFunction,
                                PotentialMinimum const& falseVacuum,
                                PotentialMinimum const& trueVacuum,
                                double const tunnelingTemperature,
                                double const actionThreshold,
                                double const requiredVacuumSeparationSquared );
  };




  // This returns either the dimensionless bounce action integrated over four
  // dimensions (for zero temperature) or the dimensionful bounce action
  // integrated over three dimensions (for non-zero temperature) for
  // tunneling from falseVacuum to trueVacuum at temperature
  // tunnelingTemperature, or an upper bound if the upper bound drops below
  // actionThreshold during the course of the calculation. The vacua are
  // assumed to already be the minima at tunnelingTemperature.
  inline double BounceAlongPathWithThreshold::BounceAction(
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                            double const tunnelingTemperature )
  {
    double actionThreshold( lnOfThermalIntegrationFactor );
    // We assume that the threshold should be the naive threshold for thermal
    // tunneling, and then take the T = 0 version if the temperature is not a
    // valid temperature (i.e. we treat T < 0 and T "= NaN" as T = 0.0).
    if( !(tunnelingTemperature > 0.0 ) )
    {
      actionThreshold = ( log( -SolitonicFactor( potentialFunction,
                                                 falseVacuum,
                                                 trueVacuum )
                                * fourVolumeOfKnownUniverseOverGevFourth
                               / log( survivalProbabilityThreshold ) ) );
    }
    return BoundedBounceAction( potentialFunction,
                                falseVacuum,
                                trueVacuum,
                                tunnelingTemperature,
                                actionThreshold,
                                ( vacuumSeparationFractionSquared
                              * falseVacuum.SquareDistanceTo( trueVacuum ) ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEALONGPATHWITHTHRESHOLD_HPP_ */
