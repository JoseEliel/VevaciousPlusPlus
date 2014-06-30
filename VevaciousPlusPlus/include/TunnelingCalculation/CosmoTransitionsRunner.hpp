/*
 * CosmoTransitionsRunner.hpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef COSMOTRANSITIONSRUNNER_HPP_
#define COSMOTRANSITIONSRUNNER_HPP_

#include "CommonIncludes.hpp"
#include "PotentialEvaluation/PotentialFunctions/PotentialFunction.hpp"
#include "PotentialEvaluation/PotentialFunctions/IWritesPythonPotential.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "BounceWithSplines.hpp"
#include "PotentialMinimization/MinuitMinimization.hpp"
#include "ThermalActionFitter.hpp"


namespace VevaciousPlusPlus
{

  class CosmoTransitionsRunner : public BounceWithSplines
  {
  public:
    CosmoTransitionsRunner( IWritesPythonPotential& pythonPotential,
                            PotentialFunction& potentialFunction,
                            TunnelingStrategy const tunnelingStrategy,
                            double const survivalProbabilityThreshold,
                            std::string const& pathToCosmotransitions,
                            unsigned int const resolutionOfDsbVacuum = 10 );
    virtual
    ~CosmoTransitionsRunner();


    // This doesn't do anything here.
    virtual void
    UpdateSelfForNewSlha( SlhaManager const& slhaManager ){}


  protected:
    static std::string pythonPotentialFilenameBase;

    IWritesPythonPotential& pythonPotential;
    std::string const pathToCosmotransitions;
    size_t const resolutionOfDsbVacuum;


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
    // use to calculate the bounce action at tunnelingTemperature.
    virtual double BounceAction( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 double const tunnelingTemperature ) const;

    // This calculates the evaporation and critical temperatures, then writes
    // and runs a Python program using the potential from
    // pythonPotentialFilename for CosmoTransitions to get an estimate of the
    // thermal dependence of the action, then uses Minuit2 to find the
    // optimal tunneling temperature, then writes and runs another Python
    // program to use CosmoTransitions to calculate the thermal action at this
    // optimal temperature.
    virtual void
    CalculateThermalTunneling( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum );
  };

} /* namespace VevaciousPlusPlus */
#endif /* COSMOTRANSITIONSRUNNER_HPP_ */
