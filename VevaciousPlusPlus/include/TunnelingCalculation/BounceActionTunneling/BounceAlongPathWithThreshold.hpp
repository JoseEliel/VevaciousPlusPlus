/*
 * BounceAlongPathWithThreshold.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEALONGPATHWITHTHRESHOLD_HPP_
#define BOUNCEALONGPATHWITHTHRESHOLD_HPP_

#include "CommonIncludes.hpp"
#include "../BounceActionTunneler.hpp"

namespace VevaciousPlusPlus
{

  class BounceAlongPathWithThreshold : public BounceActionTunneler
  {
  public:
    BounceAlongPathWithThreshold( PotentialFunction& potentialFunction,
                                  std::string const& xmlArguments );
    virtual
    ~BounceAlongPathWithThreshold();


  protected:
    BouncePathFinder* pathFinder;
    BounceActionCalculator* actionCalculator;
    double actionThreshold;
    size_t thermalIntegrationResolution;


    // This should return either the dimensionless bounce action integrated
    // over four dimensions (for zero temperature) or the dimensionful bounce
    // action integrated over three dimensions (for non-zero temperature) for
    // tunneling from falseVacuum to trueVacuum at temperature
    // tunnelingTemperature.
    virtual double BounceAction( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 double const tunnelingTemperature ) const;

    // This should set thermalSurvivalProbability and
    // dominantTemperatureInGigaElectronVolts appropriately.
    virtual void
    CalculateThermalTunneling( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum );
  };




  // This returns either the dimensionless bounce action integrated over four
  // dimensions (for zero temperature) or the dimensionful bounce action
  // integrated over three dimensions (for non-zero temperature) for tunneling
  // from falseVacuum to trueVacuum at temperature tunnelingTemperature, or an
  // upper bound if the upper bound drops below actionThreshold during the
  // course of the calculation. The vacua are assumed to already be the minima
  // at tunnelingTemperature.
  inline double BounceAlongPathWithThreshold::BounceAction(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                      double const tunnelingTemperature ) const
  {
    if( tunnelingTemperature <= 0.0 )
    {
      double const squareRootOfSolitonicFactor(
                potentialFunction.ScaleSquaredRelevantToTunneling( falseVacuum,
                                                                trueVacuum ) );
      actionThreshold = log( -( squareRootOfSolitonicFactor
                                * squareRootOfSolitonicFactor
                                * fourVolumeOfKnownUniverseOverGevFourth )
                              / log( survivalProbabilityThreshold ) );
    }
    pathFinder->SetVacua( falseVacuum,
                          trueVacuum );
    double bounceAction( (*actionCalculator)( pathFinder->CurrentPath() ) );
    while( ( bounceAction > actionThreshold )
           &&
           pathFinder->PathCanBeImproved() )
    {
      pathFinder->ImprovePath();
      bounceAction = (*actionCalculator)( pathFinder->CurrentPath() );
    }
    return bounceAction;
  }

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEALONGPATHWITHTHRESHOLD_HPP_ */
