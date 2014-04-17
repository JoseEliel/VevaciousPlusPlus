/*
 * CosmoTransitionsRunner.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  CosmoTransitionsRunner::CosmoTransitionsRunner(
                                 IWritesPythonPotential& potentialFunction,
                                     TunnelingStrategy const tunnelingStrategy,
                                     double const survivalProbabilityThreshold,
                                  std::string const& pathToCosmotransitions ) :
    BounceWithSplines( tunnelingStrategy,
                       survivalProbabilityThreshold ),
    pathToCosmotransitions( pathToCosmotransitions )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "CosmoTransitionsRunner::CosmoTransitionsRunner( ... )";
    std::cout << std::endl;/**/
  }

  CosmoTransitionsRunner::~CosmoTransitionsRunner()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "CosmoTransitionsRunner::~CosmoTransitionsRunner()";
    std::cout << std::endl;/**/
  }


  // This creates an entire Python program to run CosmoTransitions to
  // calculate the bounce action(s) for the tunneling strategy.
  void CosmoTransitionsRunner::CalculateTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "CosmoTransitionsRunner::CalculateTunneling( ... )";
    std::cout << std::endl;/**/

    if( ( tunnelingStrategy == NoTunneling )
        ||
        ( tunnelingStrategy == NotSet ) )
    {
      quantumSurvivalProbability = -1.0;
      quantumLifetimeInSeconds = -1.0;
      thermalSurvivalProbability = -1.0;
      dominantTemperatureInGigaElectronVolts = -1.0;
      return;
    }

    // If we should tunnel, we should write the potential in Python:



    // if quantum:
    // 1) C++: write Py potential
    // 2) Py: deform at T=0

    // if thermal:
    // 1) C++: evaporate DSB - possibly exclude based on DSB being less deep
    //                         than origin
    // 2) C++: find critical T, choose T_i vector, write vacuaPair_i vector
    // 3) C++: write Py potential
    // 4) Py: find direct S_3(T_i) vector
    // 5) C++: read S_3(T_i), fit to function, Minuit2 on function -> T_opt
    // 6) Py: deform at T_opt
  }

} /* namespace VevaciousPlusPlus */
