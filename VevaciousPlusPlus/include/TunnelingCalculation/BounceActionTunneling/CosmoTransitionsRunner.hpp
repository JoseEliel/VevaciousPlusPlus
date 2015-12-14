/*
 * CosmoTransitionsRunner.hpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef COSMOTRANSITIONSRUNNER_HPP_
#define COSMOTRANSITIONSRUNNER_HPP_

#include "TunnelingCalculation/BounceActionTunneler.hpp"
#include "TunnelingCalculation/TunnelingCalculator.hpp"
#include <string>
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include <vector>
#include <fstream>
#include <iomanip>
#include "VersionInformation.hpp"
#include <cstdlib>
#include <iostream>
#include "ThermalActionFitter.hpp"
#include "Minuit2/MnMigrad.h"
#include "PotentialMinimization/GradientBasedMinimization/MinuitPotentialMinimizer.hpp"
#include <cmath>
#include "BounceActionEvaluation/BubbleShootingOnPathInFieldSpace.hpp"
#include "BounceActionEvaluation/PathParameterization/LinearSplineThroughNodes.hpp"
#include "BounceActionEvaluation/SplinePotential.hpp"
#include "BounceActionEvaluation/BubbleProfile.hpp"
#include <limits>

namespace VevaciousPlusPlus
{

  class CosmoTransitionsRunner : public BounceActionTunneler
  {
  public:
    CosmoTransitionsRunner(
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                            double const survivalProbabilityThreshold,
                            unsigned int const temperatureAccuracy,
                            std::string const& pathToCosmotransitions,
                            unsigned int const resolutionOfDsbVacuum,
                            unsigned int const maxInnerLoops,
                            unsigned int const maxOuterLoops,
                           unsigned int const thermalStraightPathFitResolution,
                            double const vacuumSeparationFraction );
    virtual ~CosmoTransitionsRunner();


  protected:
    static std::string pythonPotentialFilenameBase;

    std::string const pathToCosmotransitions;
    unsigned int const resolutionOfDsbVacuum;
    unsigned int const maxInnerLoops;
    unsigned int const maxOuterLoops;
    unsigned int const thermalStraightPathFitResolution;


    // This creates a Python file with the potential in a form that can be used
    // by CosmoTransitions.
    virtual void
    PrepareCommonExtras( PotentialFunction const& potentialFunction )
    { potentialFunction.WriteAsPython( pythonPotentialFilenameBase + ".py" ); }

    // This returns either the dimensionless bounce action integrated over four
    // dimensions (for zero temperature) or the dimensionful bounce action
    // integrated over three dimensions (for non-zero temperature) for
    // tunneling from falseVacuum to trueVacuum at temperature
    // tunnelingTemperature. It does so by writing and running a Python program
    // using the potential from pythonPotentialFilename for CosmoTransitions to
    // use to calculate the bounce action at tunnelingTemperature. The vacua
    // are assumed to already be the minima at tunnelingTemperature.
    virtual double BounceAction( PotentialFunction const& potentialFunction,
                                 PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 double const tunnelingTemperature );

    // This calculates the evaporation and critical temperatures, then writes
    // and runs a Python program using the potential from
    // pythonPotentialFilename for CosmoTransitions to get an estimate of the
    // thermal dependence of the action, then uses Minuit2 to find the optimal
    // tunneling temperature, then writes and runs another Python program to
    // use CosmoTransitions to calculate the thermal action at this optimal
    // temperature.
    virtual void
    ContinueThermalTunneling( PotentialFunction const& potentialFunction,
                              PotentialMinimum const& falseVacuum,
                              PotentialMinimum const& trueVacuum,
                             double const potentialAtOriginAtZeroTemperature );

    // This is just to make it easy to switch back to using CosmoTransitions
    // again if somehow the internal straight path bounce action calculation
    // turns out to be deficient. It uses the given potentialFunction rather
    // than pythonPotential.
    void
    CalculateStraightPathActions( PotentialFunction const& potentialFunction,
                                  PotentialMinimum const& falseVacuum,
                                  PotentialMinimum const& trueVacuum,
                                  std::vector< double > const& fitTemperatures,
                                  std::vector< double >& straightPathActions )
    { InteralGuessFromStraightPaths( potentialFunction,
                                     falseVacuum,
                                     trueVacuum,
                                     fitTemperatures,
                                     straightPathActions ); }

    // This uses a BubbleShootingOnSpline object at different temperatures to
    // fill straightPathActions based on straight paths between the thermal
    // vacua. It uses the given potentialFunction rather than pythonPotential.
    void
    InteralGuessFromStraightPaths( PotentialFunction const& potentialFunction,
                                   PotentialMinimum const& falseVacuum,
                                   PotentialMinimum const& trueVacuum,
                                  std::vector< double > const& fitTemperatures,
                                  std::vector< double >& straightPathActions );

    // This writes a Python program using CosmoTransitions with the minimum
    // number of deformations to get a set of actions at temperatures, and then
    // reads in the file created to fill straightPathActions. The minima at
    // the various temperatures are determined with potentialFunction, while
    // the tunneling is calculated with the Python code generated by
    // pythonPotential.
    void FitFromCosmoTransitionsStraightPaths(
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                  std::vector< double > const& fitTemperatures,
                                  std::vector< double >& straightPathActions );
  };

} /* namespace VevaciousPlusPlus */
#endif /* COSMOTRANSITIONSRUNNER_HPP_ */
