/*
 * GradientFromStartingPoints.hpp
 *
 *  Created on: Jun 30, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef GRADIENTFROMSTARTINGPOINTS_HPP_
#define GRADIENTFROMSTARTINGPOINTS_HPP_

#include "PotentialMinimizer.hpp"
#include "StartingPointFinder.hpp"
#include "GradientMinimizer.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimum.hpp"
#include <vector>
#include <iostream>
#include <cmath>

namespace VevaciousPlusPlus
{
  class GradientFromStartingPoints : public PotentialMinimizer
  {
  public:
    GradientFromStartingPoints( PotentialFunction& potentialFunction,
                                std::unique_ptr<StartingPointFinder> startingPointFinder,
                                std::unique_ptr<GradientMinimizer> gradientMinimizer,
                              double const extremumSeparationThresholdFraction,
                                double const nonDsbRollingToDsbScalingFactor, bool global_Is_Panic );
    virtual ~GradientFromStartingPoints();


    // This first sets dsbVacuum from the input recorded in
    // potentialFunction.DsbFieldValues() using gradientMinimizer, then uses
    // startingPointFinder to find the starting points for gradientMinimizer,
    // then uses gradientMinimizer to minimize the potential at a temperature
    // given by minimizationTemperature, recording the found minima in
    // foundMinima. It also records the minima lower than dsbVacuum in
    // panicVacua, and of those, it sets panicVacuum to be the minimum in
    // panicVacua closest to dsbVacuum.
    virtual void FindMinima( double const minimizationTemperature = 0.0 );

    // This uses gradientMinimizer to find the minimum at temperature
    // minimizationTemperature nearest to minimumToAdjust (which is assumed to
    // be a minimum of the potential at a different temperature).
    virtual PotentialMinimum
    AdjustMinimumForTemperature( PotentialMinimum const& minimumToAdjust,
                                 double const minimizationTemperature );

    virtual void setWhichPanicVacuum( bool global_Is_Panic_setting);

  protected:
    std::unique_ptr<StartingPointFinder> startingPointFinder;
    std::unique_ptr<GradientMinimizer> gradientMinimizer;
    std::vector< std::vector< double > > startingPoints;
    double extremumSeparationThresholdFraction;
    double nonDsbRollingToDsbScalingFactor;
    bool global_Is_Panic;
    bool done_homotopy;
  };




  // This uses gradientMinimizer to find the minimum at temperature
  // minimizationTemperature nearest to minimumToAdjust (which is assumed to
  // be a minimum of the potential at a different temperature).
  inline PotentialMinimum
  GradientFromStartingPoints::AdjustMinimumForTemperature(
                                       PotentialMinimum const& minimumToAdjust,
                                         double const minimizationTemperature )
  {
    gradientMinimizer->SetTemperature( minimizationTemperature );
    return (*gradientMinimizer)( minimumToAdjust.FieldConfiguration() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* GRADIENTFROMSTARTINGPOINTS_HPP_ */
