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
#include "Minuit2/MnPrint.h"
#include "ModifiedBounceForMinuit.hpp"

namespace VevaciousPlusPlus
{

  class MinuitBounceActionMinimizer : public BounceWithSplines
  {
  public:
    MinuitBounceActionMinimizer( PotentialFunction& potentialFunction,
                                 TunnelingStrategy const tunnelingStrategy,
                                 double const survivalProbabilityThreshold,
                                 size_t const numberOfNodes,
                                 double const initialStepSize,
                                 unsigned int const minuitStrategy = 1,
                                 double const minuitTolerance = 100.0 );
    virtual
    ~MinuitBounceActionMinimizer();


    // This should perform all relevant updates for the new SLHA data except
    // for propagating the push to the set of dependent SlhaUpdatePropagators.
    virtual void
    UpdateSelfForNewSlha( SlhaManager const& slhaManager );


  protected:
    size_t const numberOfNodes;
    double const initialStepSize;
    unsigned int const minuitStrategy;
    double const minuitTolerance;


    // This returns either the dimensionless bounce action integrated over
    // four dimensions (for zero temperature) or the dimensionful bounce
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

} /* namespace VevaciousPlusPlus */
#endif /* MINUITBOUNCEACTIONMINIMIZER_HPP_ */
