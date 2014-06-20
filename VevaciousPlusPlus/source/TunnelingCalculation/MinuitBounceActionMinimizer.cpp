/*
 * MinuitBounceActionMinimizer.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  MinuitBounceActionMinimizer::MinuitBounceActionMinimizer(
                                          PotentialFunction& potentialFunction,
                                     TunnelingStrategy const tunnelingStrategy,
                                  double const survivalProbabilityThreshold ) :
    BounceWithSplines( potentialFunction,
                       tunnelingStrategy,
                       survivalProbabilityThreshold )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::MinuitBounceActionMinimizer( ... )";
    std::cout << std::endl;/**/
  }

  MinuitBounceActionMinimizer::~MinuitBounceActionMinimizer()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::~MinuitBounceActionMinimizer()";
    std::cout << std::endl;/**/
  }


  // This doesn't do anything here.
  void MinuitBounceActionMinimizer::UpdateSelfForNewSlha(
                                               SlhaManager const& slhaManager )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::UpdateSelfForNewSlha(...)";
    std::cout << std::endl;/**/
  }

  // This should set quantumSurvivalProbability and quantumLifetimeInSeconds
  // appropriately.
  void MinuitBounceActionMinimizer::CalculateQuantumTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::CalculateQuantumTunneling(...)";
    std::cout << std::endl;/**/

    criticalRatherThanEvaporation = false;
    evaporationMinimum = falseVacuum;
    //double evaporationTemperature( CriticalOrEvaporationTemperature(
    //          potentialFunction( potentialFunction.FieldValuesOrigin() ) ) );
    ModifiedBounceForMinuit modifiedBounceForMinuit( potentialFunction,
                                                     3,
                                                     16,
                                                     falseVacuum,
                                                     trueVacuum,
                                                     0.0,
                                                     32 );
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "Remember to write in taking shoot attempts, path resolution, etc.,"
    << " from input XML.";
    std::cout << std::endl;/**/
    std::vector< double > pathParameterization;
    std::vector< double > initialStepSizes;
    modifiedBounceForMinuit.SetUpStraightPathForMinuit( pathParameterization,
                                                        initialStepSizes,
                                                        0.0,
                                                        0.5 );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "straight path bounce = "
    << modifiedBounceForMinuit( pathParameterization );
    std::cout << std::endl;
    if( !(initialStepSizes.empty()) )
    {
      initialStepSizes.back() = 0.0;
      std::cout
      << "test path bounce = " << modifiedBounceForMinuit( initialStepSizes );
    }
    std::cout << std::endl;/**/

    // Should now do migrad...
  }

  // This should set thermalSurvivalProbability and
  // dominantTemperatureInGigaElectronVolts appropriately.
  void MinuitBounceActionMinimizer::CalculateThermalTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::CalculateThermalTunneling(...)";
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
