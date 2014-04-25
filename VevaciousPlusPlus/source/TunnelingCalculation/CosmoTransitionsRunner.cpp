/*
 * CosmoTransitionsRunner.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{
  std::string
  CosmoTransitionsRunner::pythonPotentialFilenameBase( "VevaciousPotential");

  CosmoTransitionsRunner::CosmoTransitionsRunner(
                                       IWritesPythonPotential& pythonPotential,
                                          PotentialFunction& potentialFunction,
                                     TunnelingStrategy const tunnelingStrategy,
                                     double const survivalProbabilityThreshold,
                                  std::string const& pathToCosmotransitions ) :
    BounceWithSplines( potentialFunction,
                       tunnelingStrategy,
                       survivalProbabilityThreshold ),
    pythonPotential( pythonPotential ),
    pathToCosmotransitions( pathToCosmotransitions )
  {
    // This constructor is just an initialization list.
  }

  CosmoTransitionsRunner::~CosmoTransitionsRunner()
  {
    // This does nothing.
  }



  // This writes and runs a Python program using the potential from
  // pythonPotentialFilename for CosmoTransitions.
  void CosmoTransitionsRunner::CalculateQuantumTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    std::string pythonResultFilename( "VevaciousCosmoTransitionsResult.txt" );
    std::string pythonMainFilename( "VevaciousCosmoTransitionsRunner.py" );
    std::ofstream pythonFile;
    pythonFile.open( pythonMainFilename.c_str() );
    pythonFile << std::setprecision( 12 );
    pythonFile << "# Automatically generated file! Modify at your own PERIL!\n"
    "# This file was created by Vevacious version "
    << VevaciousPlusPlus::versionString << "\n"
    "from __future__ import division\n"
    "import sys\n"
    "import math\n"
    "import numpy\n"
    "import " << pythonPotentialFilenameBase << " as VPD\n"
    "\n"
    "pathToCosmotransitions = \"" << pathToCosmotransitions << "\"\n"
    "sys.path.append( pathToCosmotransitions )\n"
    "import pathDeformation as CTPD\n"
    "trueAndFalseVacua = [ [ ";
    for( std::vector< double >::const_iterator
         fieldValue( trueVacuum.FieldConfiguration().begin() );
         fieldValue < trueVacuum.FieldConfiguration().end();
         ++fieldValue )
    {
      if( fieldValue != trueVacuum.FieldConfiguration().begin() )
      {
        pythonFile << ", ";
      }
      pythonFile << *fieldValue;
    }
    pythonFile << " ],\n"
    "[ ";
    for( std::vector< double >::const_iterator
         fieldValue( falseVacuum.FieldConfiguration().begin() );
         fieldValue < falseVacuum.FieldConfiguration().end();
         ++fieldValue )
    {
      if( fieldValue != falseVacuum.FieldConfiguration().begin() )
      {
        pythonFile << ", ";
      }
      pythonFile << *fieldValue;
    }
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "trueVacuum.SquareDistanceTo( falseVacuum ) = "
    << trueVacuum.SquareDistanceTo( falseVacuum )
    << std::endl
    << "trueVacuum.LengthSquared() = "
    << trueVacuum.LengthSquared()
    << std::endl
    << "sqrt( trueVacuum.SquareDistanceTo( falseVacuum )"
    << " / trueVacuum.LengthSquared() ) = "
    << sqrt( trueVacuum.SquareDistanceTo( falseVacuum )
             / trueVacuum.LengthSquared() );
    std::cout << std::endl;/**/

    pythonFile << " ] ]\n"
    "tunnelPathPoints = "
    << (int)( 10.0 * sqrt( trueVacuum.SquareDistanceTo( falseVacuum )
                           / trueVacuum.LengthSquared() ) ) << "\n"
    "VPD.UnderlyingPotential = VPD.JustLoopCorrectedPotential\n"
    "tunnelingResult = CTPD.fullTunneling( path_pts = trueAndFalseVacua,\n"
    "                                  V = VPD.PotentialForCosmotransitions,\n"
    "                                  dV = VPD.GradientForCosmotransitions,\n"
    "                             tunneling_init_params = dict( alpha = 3 ),\n"
    "     tunneling_findProfile_params = dict( npoints = tunnelPathPoints ),\n"
    "                      deformation_deform_params = dict( maxiter = 10 ),\n"
    "                                      maxiter = 10 )\n"
    "tunnelingCalculator.run( maxiter2 = 3 )\n"
    "# maxiter2 is a keyword argument which gets passed along to the\n"
    "# CosmoTransitions fullTunneling object when it calls self.run(), which\n"
    "# controls how many iterations of the outer loop of the deformation\n"
    "# happen.\n"
    "outputFile = open( \"" << pythonResultFilename << "\", \"w\" )\n"
    "outputFile.write( str( tunnelingCalculator.findAction() ) )\n"
    "outputFile.close()\n"
    "\n"
    "# End of automatically generated file!\n";


    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "Need to do something about barriers smaller than tunnel path"
    << " resolution!";
    std::cout << std::endl;/**/

    pythonFile.close();

    BOL::UsefulStuff::runSystemCommand( "python " + pythonMainFilename );

    double quantumAction( -1.0 );
    std::ifstream resultStream;
    resultStream.open( pythonResultFilename.c_str() );
    resultStream >> quantumAction;
    resultStream.close();

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "quantumAction = " << quantumAction;
    std::cout << std::endl;/**/

    if( quantumAction >= maximumPowerOfNaturalExponent )
    {
      std::cout
      << std::endl
      << "Warning! The calculated bounce action was so large that"
      << " exponentiating it would result in an overflow error, so capping it"
      << " at " << maximumPowerOfNaturalExponent;
      std::cout << std::endl;
      quantumAction = maximumPowerOfNaturalExponent;
    }
    else if( quantumAction <= -maximumPowerOfNaturalExponent )
    {
      std::cout
      << std::endl
      << "Warning! The calculated bounce action was so large and negative that"
      << " exponentiating it would result in an overflow error, so capping it"
      << " at " << -maximumPowerOfNaturalExponent;
      std::cout << std::endl;
      quantumAction = -maximumPowerOfNaturalExponent;
    }
    quantumLifetimeInSeconds
    = ( ( exp( 0.25 * quantumAction )
          * hBarInGigaElectronVoltSeconds )
        / sqrt( potentialFunction.ScaleSquaredRelevantToTunneling( falseVacuum,
                                                              trueVacuum ) ) );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "quantumLifetimeInSeconds = " << quantumLifetimeInSeconds;
    std::cout << std::endl;/**/

    double
    survivalExponent( ageOfKnownUniverseInSeconds / quantumLifetimeInSeconds );
    if( survivalExponent >= maximumPowerOfNaturalExponent )
    {
      std::cout
      << std::endl
      << "Warning! The calculated decay width was so large that"
      << " exponentiating it would result in an overflow error, so setting the"
      << " survival probability to zero.";
      quantumSurvivalProbability = 0.0;
    }
    else
    {
      quantumSurvivalProbability = exp( -survivalExponent );
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "quantumSurvivalProbability = " << quantumSurvivalProbability;
    std::cout << std::endl;/**/
  }

  // This calculates the evaporation and critical temperatures, then writes
  // and runs a Python program using the potential from
  // pythonPotentialFilename for CosmoTransitions to get an estimate of the
  // thermal dependence of the action, then uses Minuit2 to find the
  // optimal tunneling temperature, then writes and runs another Python
  // program to use CosmoTransitions to calculate the thermal action at this
  // optimal temperature.
  void CosmoTransitionsRunner::CalculateThermalTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // 1) C++: evaporate DSB - possibly exclude based on DSB being less deep
    //                         than origin
    // 2) C++: find critical T, choose T_i vector, write vacuaPair_i vector
    // 3) C++: write Py potential
    // 4) Py: find direct S_3(T_i) vector
    // 5) C++: read S_3(T_i), fit to function, Minuit2 on function -> T_opt
    // 6) Py: deform at T_opt
  }

} /* namespace VevaciousPlusPlus */
