/*
 * HomotopyContinuationAndMinuit.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONANDMINUIT_HPP_
#define HOMOTOPYCONTINUATIONANDMINUIT_HPP_

#include "../CommonIncludes.hpp"
#include "HomotopyContinuationAndGradient.hpp"
#include "../PotentialEvaluation.hpp"
#include "MinuitMinimization.hpp"

namespace VevaciousPlusPlus
{

  class HomotopyContinuationAndMinuit : public HomotopyContinuationAndGradient
  {
  public:
    HomotopyContinuationAndMinuit( PotentialFunction const& potentialFunction,
                        HomotopyContinuationSolver& homotopyContinuationSolver,
                            double const extremumSeparationThresholdFraction );
    virtual
    ~HomotopyContinuationAndMinuit();


    // This finds all the extrema of polynomialPotential with
    // homotopyContinuationSolver, then uses them as starting points for
    // Minuit2 to minimize potentialForMinuit, evaluated at a temperature given
    // by temperatureInGev, and records the found minima in foundMinima. It
    // also sets dsbVacuum (using polynomialPotential.DsbFieldValues() as a
    // starting point for Minuit2), and records the minima lower than dsbVacuum
    // in panicVacua, and of those, it sets panicVacuum to be the minimum in
    // panicVacua closest to dsbVacuum.
    virtual void FindMinima( double const temperatureInGev = 0.0 );

    // This should find the minimum at temperature temperatureInGev nearest to
    // minimumToAdjust (which is assumed to be a minimum of the potential at a
    // different temperature).
    virtual PotentialMinimum
    AdjustMinimumForTemperature( PotentialMinimum const& minimumToAdjust,
                                 double const temperatureInGev );


  protected:
    PotentialForMinuit potentialForMinuit;
    MinuitManager minuitManager;
    double const extremumSeparationThresholdFraction;

    // This uses Minuit2 to minimize potentialForMinuit starting from the
    // values in purelyRealSolutionSets.
    void RollAndSortExtrema();
  };




  // This finds all the extrema of polynomialPotential with
  // homotopyContinuationSolver, then uses them as starting points for
  // Minuit2 to minimize potentialForMinuit, evaluated at a temperature given
  // by temperatureInGev, and records the found minima in foundMinima. It
  // also sets dsbVacuum (using polynomialPotential.DsbFieldValues() as a
  // starting point for Minuit2), and records the minima lower than dsbVacuum
  // in panicVacua, and of those, it sets panicVacuum to be the minimum in
  // panicVacua closest to dsbVacuum.
  inline void
  HomotopyContinuationAndMinuit::FindMinima( double const temperatureInGev )
  {
    homotopyContinuationSolver.FindTreeLevelExtrema( purelyRealSolutionSets );
    potentialForMinuit.SetTemperature( temperatureInGev );
    dsbVacuum = minuitManager( potentialFunction.DsbFieldValues() );
    RollAndSortExtrema();
  }

  // This should find the minimum at temperature temperatureInGev nearest to
  // minimumToAdjust (which is assumed to be a minimum of the potential at a
  // different temperature).
  inline PotentialMinimum
  HomotopyContinuationAndMinuit::AdjustMinimumForTemperature(
                                       PotentialMinimum const& minimumToAdjust,
                                                double const temperatureInGev )
  {
    potentialForMinuit.SetTemperature( temperatureInGev );
    return
    PotentialMinimum( minuitManager( minimumToAdjust.FieldConfiguration() ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONANDMINUIT_HPP_ */
