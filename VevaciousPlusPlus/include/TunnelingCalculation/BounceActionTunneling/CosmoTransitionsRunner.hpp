/*
 * CosmoTransitionsRunner.hpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef COSMOTRANSITIONSRUNNER_HPP_
#define COSMOTRANSITIONSRUNNER_HPP_

#include "CommonIncludes.hpp"
#include "limits"
#include "VersionInformation.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialEvaluation/PotentialFunctions/IWritesPythonPotential.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "../BounceActionTunneler.hpp"
#include "ThermalActionFitter.hpp"
#include "BounceActionEvaluation/PathParameterization/LinearSplineThroughNodes.hpp"
#include "BounceActionEvaluation/BubbleShootingOnSpline.hpp"
#include "BounceActionEvaluation/BubbleProfile.hpp"

namespace VevaciousPlusPlus
{

  class CosmoTransitionsRunner : public BounceActionTunneler
  {
  public:
    CosmoTransitionsRunner( IWritesPythonPotential& pythonPotential,
                            PotentialFunction& potentialFunction,
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                            double const survivalProbabilityThreshold,
                            size_t const temperatureAccuracy,
                            std::string const& pathToCosmotransitions,
                            size_t const resolutionOfDsbVacuum,
                            size_t const maxInnerLoops,
                            size_t const maxOuterLoops,
                            size_t const thermalStraightPathFitResolution );
    virtual ~CosmoTransitionsRunner();


  protected:
    static std::string pythonPotentialFilenameBase;

    IWritesPythonPotential& pythonPotential;
    std::string const pathToCosmotransitions;
    size_t const resolutionOfDsbVacuum;
    size_t const maxInnerLoops;
    size_t const maxOuterLoops;
    size_t const thermalStraightPathFitResolution;


    // This creates a Python file with the potential in a form that can be used
    // by CosmoTransitions.
    virtual void PrepareCommonExtras()
    { pythonPotential.WriteAsPython( pythonPotentialFilenameBase + ".py" ); }

    // This returns either the dimensionless bounce action integrated over four
    // dimensions (for zero temperature) or the dimensionful bounce action
    // integrated over three dimensions (for non-zero temperature) for
    // tunneling from falseVacuum to trueVacuum at temperature
    // tunnelingTemperature. It does so by writing and running a Python program
    // using the potential from pythonPotentialFilename for CosmoTransitions to
    // use to calculate the bounce action at tunnelingTemperature. The vacua
    // are assumed to already be the minima at tunnelingTemperature.
    virtual double BounceAction( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 double const tunnelingTemperature );

    // This calculates the evaporation and critical temperatures, then writes
    // and runs a Python program using the potential from
    // pythonPotentialFilename for CosmoTransitions to get an estimate of the
    // thermal dependence of the action, then uses Minuit2 to find the
    // optimal tunneling temperature, then writes and runs another Python
    // program to use CosmoTransitions to calculate the thermal action at this
    // optimal temperature.
    virtual void ContinueThermalTunneling( PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum,
                             double const potentialAtOriginAtZeroTemperature );


    // This is just to make it easy to switch back to using CosmoTransitions
    // again if somehow the internal straight path bounce action calculation
    // turns out to be deficient.
    void
    CalculateStraightPathActions( PotentialMinimum const& falseVacuum,
                                  PotentialMinimum const& trueVacuum,
                                  std::vector< double > const& fitTemperatures,
                                  std::vector< double >& straightPathActions )
    { InteralGuessFromStraightPaths( falseVacuum,
                                     trueVacuum,
                                     fitTemperatures,
                                     straightPathActions ); }

    // This uses a BubbleShootingOnSpline object at different temperatures to
    // fill straightPathActions based on straight paths between the thermal
    // vacua.
    void InteralGuessFromStraightPaths( PotentialMinimum const& falseVacuum,
                                        PotentialMinimum const& trueVacuum,
                                  std::vector< double > const& fitTemperatures,
                                  std::vector< double >& straightPathActions );

    // This writes a Python programme using CosmoTransitions with the minimum
    // number of deformations to get a set of actions at temperatures, and then
    // reads in the file created to fill straightPathActions.
    void
    FitFromCosmoTransitionsStraightPaths( PotentialMinimum const& falseVacuum,
                                          PotentialMinimum const& trueVacuum,
                                  std::vector< double > const& fitTemperatures,
                                  std::vector< double >& straightPathActions );
  };

} /* namespace VevaciousPlusPlus */
#endif /* COSMOTRANSITIONSRUNNER_HPP_ */
