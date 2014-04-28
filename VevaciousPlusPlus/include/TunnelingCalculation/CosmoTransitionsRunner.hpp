/*
 * CosmoTransitionsRunner.hpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef COSMOTRANSITIONSRUNNER_HPP_
#define COSMOTRANSITIONSRUNNER_HPP_

#include "../CommonIncludes.hpp"
#include "../PotentialEvaluation.hpp"
#include "../PotentialMinimization.hpp"
#include "BounceWithSplines.hpp"
#include "../PotentialMinimization/MinuitMinimization.hpp"
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
                            std::string const& pathToCosmotransitions );
    virtual
    ~CosmoTransitionsRunner();


    // This doesn't do anything here.
    virtual void
    UpdateSelfForNewSlha( SlhaManager const& slhaManager ){}


  protected:
    static std::string pythonPotentialFilenameBase;
    IWritesPythonPotential& pythonPotential;
    std::string const pathToCosmotransitions;

    // This creates a Python file with the potential in a form that can be used
    // by CosmoTransitions.
    virtual void PrepareCommonExtras()
    { pythonPotential.WriteAsPython( pythonPotentialFilenameBase + ".py" ); }

    // This writes and runs a Python program using the potential from
    // pythonPotentialFilename for CosmoTransitions.
    virtual void
    CalculateQuantumTunneling( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum );

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

    double const DeformedPathAction( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum,
                                     double const tunnelingTemperature ) const;

    // This calculates the temperature at which either tunneling from
    // givenVacuum to the field origin becomes impossible if
    // criticalRatherThanEvaporation is true or the temperature at which
    // givenVacuum evaporates if false.
    double const CriticalTemperature( PotentialMinimum const& givenVacuum,
                                      bool const criticalRatherThanEvaporation,
                                      double const potentialAtOrigin ) const;

    // This returns true if the temperature is above the temperature that
    // is calculated by CriticalTemperature above.
    bool const BelowCritical( bool const criticalRatherThanEvaporation,
                              PotentialMinimum const& thermalMinimum,
                              double const thresholdSeparationSquared,
                              double const temperatureGuess ) const;
  };



  // This returns true if the temperature is above the temperature that
  // is calculated by CriticalTemperature above.
  inline bool const CosmoTransitionsRunner::BelowCritical(
                                      bool const criticalRatherThanEvaporation,
                                        PotentialMinimum const& thermalMinimum,
                                       double const thresholdSeparationSquared,
                                          double const temperatureGuess ) const
  {
    if( criticalRatherThanEvaporation )
    {
      return ( potentialFunction( thermalMinimum.FieldConfiguration(),
                                  temperatureGuess )
               < potentialFunction( potentialFunction.FieldValuesOrigin(),
                                    temperatureGuess ) );
    }
    return ( thermalMinimum.LengthSquared() > thresholdSeparationSquared );
  }

} /* namespace VevaciousPlusPlus */
#endif /* COSMOTRANSITIONSRUNNER_HPP_ */
