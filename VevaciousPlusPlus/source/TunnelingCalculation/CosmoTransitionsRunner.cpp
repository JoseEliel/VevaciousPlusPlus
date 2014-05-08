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
                                     std::string const& pathToCosmotransitions,
                                  unsigned int const resolutionOfDsbVacuum ) :
    BounceWithSplines( potentialFunction,
                       tunnelingStrategy,
                       survivalProbabilityThreshold ),
    pythonPotential( pythonPotential ),
    pathToCosmotransitions( pathToCosmotransitions ),
    resolutionOfDsbVacuum( resolutionOfDsbVacuum )
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
    double quantumAction( DeformedPathAction( falseVacuum,
                                              trueVacuum,
                                              0.0 ) );
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
    // First we find the temperature at which the DSB vacuum evaporates, and
    // possibly exclude the parameter point based on DSB being less deep than
    // origin.
    double const potentialAtOrigin( potentialFunction(
                                     potentialFunction.FieldValuesOrigin() ) );
    if( potentialFunction( falseVacuum.FieldConfiguration() )
        > potentialAtOrigin )
    {
      std::cout
      << std::endl
      << "DSB vacuum has higher energy density than vacuum with no non-zero"
      << " VEVs! Assuming that it is implausible that the Universe cooled into"
      << " this false vacuum from the symmetric phase, and so setting survival"
      << " probability to 0.";
      std::cout << std::endl;
      dominantTemperatureInGigaElectronVolts = 0.0;
      thermalSurvivalProbability = 0.0;
      return;
    }

    criticalRatherThanEvaporation = false;
    thresholdSeparationSquared = ( 0.05 * 0.05 * falseVacuum.LengthSquared() );
    evaporationMinimum = falseVacuum;
    double const
    falseEvaporationTemperature( CriticalOrEvaporationTemperature(
                                                         potentialAtOrigin ) );

    // Next we find a set of true and false vacua at a series of temperatures
    // so that we can fit the temperature dependence of the thermal action.

    criticalRatherThanEvaporation = true;
    criticalMinimum = trueVacuum;
    // We need the temperature where tunneling from the true vacuum to the
    // field origin becomes impossible.
    double const
    criticalTunnelingTemperature( CriticalOrEvaporationTemperature(
                                                         potentialAtOrigin ) );
    // We are going to fit a function that diverges at
    // criticalTunnelingTemperature, so we shouldn't pick a point too close to
    // close to the divergence.
    double const highestFitTemperature( 0.9 * criticalTunnelingTemperature );

    unsigned int const nodesForFit( 5 );
    std::vector< std::vector< double > > fitFalseVacua( nodesForFit );
    std::vector< std::vector< double > > fitTrueVacua( nodesForFit );
    std::vector< double > fitTemperatures( nodesForFit );

    // Ideally the DSB vacuum evaporates at a much lower temperature than the
    // critical tunneling temperature for the true vacuum, so we can restrict
    // our fit to temperatures above the evaporation temperature where there is
    // a kink in the thermal action due to the rapid acceleration of the false
    // vacuum end of the tunneling path followed by a complete halt, with
    // respect to increasing temperature.
    if( criticalTunnelingTemperature > ( 2.0 * falseEvaporationTemperature ) )
    {
      double const lowestFitTemperature( 1.1 * falseEvaporationTemperature );
      double const
      stepTemperature( ( highestFitTemperature - lowestFitTemperature )
                       / (double)( nodesForFit - 1 ) );
      double currentTemperature( lowestFitTemperature );
      for( unsigned int whichNode( 0 );
           whichNode < nodesForFit;
           ++whichNode )
      {
        fitTemperatures[ whichNode ] = currentTemperature;
        potentialForMinuit.SetTemperature( currentTemperature );
        fitFalseVacua[ whichNode ] = potentialFunction.DsbFieldValues();
        fitTrueVacua[ whichNode ] = thermalMinimumMinuit(
                            trueVacuum.FieldConfiguration() ).VariableValues();
        currentTemperature += stepTemperature;
      }
    }
    else
    {
      // If we have to fit below the kink, we have to do so...
      double const lowestFitTemperature( 1.0 );
      double const
      stepTemperature( ( highestFitTemperature - lowestFitTemperature )
                       / (double)( nodesForFit - 1 ) );
      double currentTemperature( lowestFitTemperature );
      for( unsigned int whichNode( 0 );
           whichNode < nodesForFit;
           ++whichNode )
      {
        fitTemperatures[ whichNode ] = currentTemperature;
        potentialForMinuit.SetTemperature( currentTemperature );
        fitFalseVacua[ whichNode ] = thermalMinimumMinuit(
                           falseVacuum.FieldConfiguration() ).VariableValues();
        fitTrueVacua[ whichNode ] = thermalMinimumMinuit(
                            trueVacuum.FieldConfiguration() ).VariableValues();
        currentTemperature += stepTemperature;
      }
    }
    std::string const
    pythonResultFilename( "VevaciousCosmoTransitionsThermalFitResult.txt" );
    std::string const
    pythonMainFilename( "VevaciousCosmoTransitionsThermalFitter.py" );
    std::string systemCommand( "rm " );
    systemCommand.append( pythonMainFilename );
    systemCommand.append( "c" );
    BOL::UsefulStuff::runSystemCommand( systemCommand );
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
    "ctVersionString = getattr( CTPD, \"__version__\", \"1\" )\n"
    "ctMajorVersion = int( ctVersionString.split( \'.\' )[ 0 ] )\n"
    "tunnelPathPoints = "
    << TunnelPathResolution( falseVacuum,
                             trueVacuum,
                             resolutionOfDsbVacuum ) << "\n"
    "\n"
    "# CosmoTransitions nests a loop within a loop to do its deformations:\n"
    "innerLoopMaxDeformations = 0\n"
    "outerLoopMaxDeformations = 1\n"
    "VPD.UnderlyingPotential = VPD.LoopAndThermallyCorrectedPotential\n"
    "tunnelingSymmetryDimensionMinusOne = 2\n"
    "\n"
    "resultString = \"error\"\n"
    "# Each tuple of fitTuples is\n"
    "# [ temperature, true vacuum array, false vacuum array].\n"
    "fitTuples = [ ";
    for( unsigned int whichNode( 0 );
         whichNode < nodesForFit;
         ++whichNode )
    {
      if( whichNode > 0 )
      {
        pythonFile << ",\n ";
      }
      pythonFile << "[ " << fitTemperatures[ whichNode ] << ", [";
      for( std::vector< double >::const_iterator
           fieldValue( fitTrueVacua[ whichNode ].begin() );
           fieldValue < fitTrueVacua[ whichNode ].end();
           ++fieldValue )
      {
        if( fieldValue != fitTrueVacua[ whichNode ].begin() )
        {
          pythonFile << ", ";
        }
        pythonFile << *fieldValue;
      }
      pythonFile << "], [";
      for( std::vector< double >::const_iterator
           fieldValue( fitFalseVacua[ whichNode ].begin() );
           fieldValue < fitFalseVacua[ whichNode ].end();
           ++fieldValue )
      {
        if( fieldValue != fitFalseVacua[ whichNode ].begin() )
        {
          pythonFile << ", ";
        }
        pythonFile << *fieldValue;
      }
      pythonFile << " ] ]";
    }
    pythonFile << " ]\n"
    "outputFile = open( \"" << pythonResultFilename << "\", \"w\" )\n"
    "for fitTuple in fitTuples:\n"
    "    VPD.invTSq = ( 1.0 / ( (fitTuple[ 0 ])**2 ) )\n"
    "    trueAndFalseVacua = [ fitTuple[ 1 ], fitTuple[ 2 ] ]\n"
    "    print( \"Calculating the thermal action along a direct path from \"\n"
    "           + str( fitTuple[ 2 ] )\n"
    "           + \" to \"\n"
    "           + str( fitTuple[ 1 ] )\n"
    "           + \" at temperature \"\n"
    "           + str( fitTuple[ 0 ] )\n"
    "           + \" GeV.\" )\n"
    "    if ( ctMajorVersion is 2 ):\n"
    "        tunnelingResult = CTPD.fullTunneling(\n"
    "                                         path_pts = trueAndFalseVacua,\n"
    "                                  V = VPD.PotentialForCosmotransitions,\n"
    "                                  dV = VPD.GradientForCosmotransitions,\n"
    "                                          tunneling_init_params = dict(\n"
    "                                              alpha = 2 ),\n"
    "                                   tunneling_findProfile_params = dict(\n"
    "                                          npoints = tunnelPathPoints ),\n"
    "                                      deformation_deform_params = dict(\n"
    "                                  maxiter = innerLoopMaxDeformations ),\n"
    "                                   maxiter = outerLoopMaxDeformations )\n"
    "        outputFile.write( \" \" + str( tunnelingResult.action ) )\n"
    "    elif ( ctMajorVersion is 1 ):\n"
    "        tunnelingCalculator = CTPD.fullTunneling(\n"
    "                                               phi = trueAndFalseVacua,\n"
    "                                  V = VPD.PotentialForCosmotransitions,\n"
    "                                  dV = VPD.GradientForCosmotransitions,\n"
    "                                                  alpha = 2,\n"
    "                                           npoints = tunnelPathPoints )\n"
    "        tunnelingCalculator.tunnel1D()\n"
    "        outputFile.write( \" \"\n"
    "                          + str( tunnelingCalculator.findAction() ) )\n"
    "\n"
    "outputFile.close()\n"
    "\n"
    "# End of automatically generated file!\n";
    pythonFile.close();
    systemCommand.assign( "python " );
    systemCommand.append( pythonMainFilename );
    std::cout
    << std::endl
    << "About to run custom Python program calling CosmoTransitions!"
    << std::endl << "Unfortunately it is likely to take quite some time (at"
    << " least 10 minutes for 4 fields at 1-loop order, probably at least an"
    << " hour for 6 fields) and the output to the terminal can lag a lot (it"
    << " might only show up after the Python has finished even)."
    << std::endl << "Calling system( \" " << systemCommand << " \")..."
    << std::endl << "-----------------";
    std::cout << std::endl;
    BOL::UsefulStuff::runSystemCommand( systemCommand );
    std::cout
    << std::endl << "-----------------" << std::endl
    << "Parsing output from " << pythonMainFilename << ".";
    std::cout << std::endl;

    std::vector< double > fittedActions( nodesForFit );
    std::ifstream resultStream;
    resultStream.open( pythonResultFilename.c_str() );
    for( unsigned int whichNode( 0 );
         whichNode < nodesForFit;
         ++whichNode )
    {
      resultStream >> fittedActions[ whichNode ];
    }
    resultStream.close();

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Fitted actions as [ T, {DSB}, {panic}, S ]:" << std::endl;
    for( unsigned int whichNode( 0 );
         whichNode < nodesForFit;
         ++whichNode )
    {
      std::cout << "[ " << fitTemperatures[ whichNode ] << ", { ";
      for( std::vector< double >::const_iterator
           fieldValue( fitFalseVacua[ whichNode ].begin() );
           fieldValue < fitFalseVacua[ whichNode ].end();
           ++fieldValue )
      {
        if( fieldValue != fitFalseVacua[ whichNode ].begin() )
        {
          std::cout << ", ";
        }
        std::cout << *fieldValue;
      }
      std::cout << "}, {";
      for( std::vector< double >::const_iterator
           fieldValue( fitTrueVacua[ whichNode ].begin() );
           fieldValue < fitTrueVacua[ whichNode ].end();
           ++fieldValue )
      {
        if( fieldValue != fitTrueVacua[ whichNode ].begin() )
        {
          std::cout << ", ";
        }
        std::cout << *fieldValue;
      }
      std::cout << "}, " << fittedActions[ whichNode ] << " ]" << std::endl;
      std::cout << std::endl;
    }
    std::cout << std::endl;/**/

    ThermalActionFitter thermalActionFitter( fitTemperatures,
                                             fittedActions,
                                             criticalTunnelingTemperature );
    MinuitManager directPathMinimizer( thermalActionFitter );
    MinuitMinimum
    directPathMinimum( directPathMinimizer( std::vector< double >( 1,
                                  ( 0.5 * criticalTunnelingTemperature ) ) ) );
    dominantTemperatureInGigaElectronVolts
    = directPathMinimum.VariableValues()[ 0 ];

    std::cout
    << std::endl
    << "Optimal temperature for tunneling estimated to be "
    << dominantTemperatureInGigaElectronVolts << " GeV.";
    std::cout << std::endl;

    // Finally we allow CosmoTransitions to calculate the action at our best
    // guess of the optimal tunneling temperature with full path deformation.
    double thermalAction( DeformedPathAction( falseVacuum,
                                              trueVacuum,
                                    dominantTemperatureInGigaElectronVolts ) );
    if( thermalAction >= ( maximumPowerOfNaturalExponent
                           * dominantTemperatureInGigaElectronVolts ) )
    {
      std::cout
      << std::endl
      << "Warning! The calculated bounce action was so large that"
      << " exponentiating it would result in an overflow error, so capping it"
      << " at " << maximumPowerOfNaturalExponent;
      std::cout << std::endl;
      thermalAction = ( maximumPowerOfNaturalExponent
                        * dominantTemperatureInGigaElectronVolts );
    }
    else if( thermalAction <= -( maximumPowerOfNaturalExponent
                                 * dominantTemperatureInGigaElectronVolts ) )
    {
      std::cout
      << std::endl
      << "Warning! The calculated bounce action was so large and negative that"
      << " exponentiating it would result in an overflow error, so capping it"
      << " at " << -maximumPowerOfNaturalExponent;
      std::cout << std::endl;
      thermalAction = -( maximumPowerOfNaturalExponent
                         * dominantTemperatureInGigaElectronVolts );
    }

    double survivalExponent( lnOfThermalIntegrationFactor
                             - ( thermalAction
                                 / dominantTemperatureInGigaElectronVolts ) );
    if( survivalExponent <= -maximumPowerOfNaturalExponent )
    {
      std::cout
      << std::endl
      << "Warning! The calculated bounce action was so large and negative that"
      << " exponentiating it would result in an overflow error, so setting the"
      << " survival probability to 1.";
      std::cout << std::endl;
      thermalSurvivalProbability = 1.0;
      return;
    }
    // We don't need to consider
    // survivalExponent > maximumPowerOfNaturalExponent as thermalAction and
    // dominantTemperatureInGigaElectronVolts should both be positive, and
    // lnOfThermalIntegrationFactor < maximumPowerOfNaturalExponent.
    double const integratedDecayWidth( exp( thermalAction ) / thermalAction );
    if( integratedDecayWidth > maximumPowerOfNaturalExponent )
    {
      std::cout
      << std::endl
      << "Warning! The calculated integrated thermal decay width was so"
      << " large that exponentiating it would result in an overflow error,"
      << " so setting the survival probability to 0.";
      std::cout << std::endl;
      thermalSurvivalProbability = 0.0;
      return;
    }
    // We don't need to consider
    // integratedDecayWidth < -maximumPowerOfNaturalExponent as
    // exp( survivalExponent ) must be positive, and thermalAction should
    // be positive.
    thermalSurvivalProbability = exp( -integratedDecayWidth );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "thermalSurvivalProbability = " << thermalSurvivalProbability;
    std::cout << std::endl;/**/
  }

  // This writes and runs a Python program using the potential from
  // pythonPotentialFilename for CosmoTransitions.
  double const CosmoTransitionsRunner::DeformedPathAction(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                      double const tunnelingTemperature ) const
  {
    std::string const
    pythonResultFilename( "VevaciousCosmoTransitionsResult.txt" );
    std::string const
    pythonMainFilename( "VevaciousCosmoTransitionsRunner.py" );
    std::string systemCommand( "rm " );
    systemCommand.append( pythonMainFilename );
    systemCommand.append( "c" );
    BOL::UsefulStuff::runSystemCommand( systemCommand );
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
    "ctVersionString = getattr( CTPD, \"__version__\", \"1\" )\n"
    "ctMajorVersion = int( ctVersionString.split( \'.\' )[ 0 ] )\n"
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
    pythonFile << " ] ]\n"
    "tunnelPathPoints = " << TunnelPathResolution( falseVacuum,
                                                   trueVacuum,
                                                resolutionOfDsbVacuum ) << "\n"
    "\n"
    "# CosmoTransitions nests a loop within a loop to do its deformations:\n"
    "innerLoopMaxDeformations = 10\n"
    "outerLoopMaxDeformations = 10\n";
    int tunnelingSymmetryDimensionMinusOne( 3 );
    std::string underlyingPotential( "JustLoopCorrectedPotential" );
    if( tunnelingTemperature > 0.0 )
    {
      pythonFile << "VPD.invTSq = "
      << ( 1.0 / ( tunnelingTemperature * tunnelingTemperature ) ) << "\n";
      tunnelingSymmetryDimensionMinusOne = 2;
      underlyingPotential.assign( "LoopAndThermallyCorrectedPotential" );
    }
    pythonFile << "tunnelingSymmetryDimensionMinusOne = "
                             << tunnelingSymmetryDimensionMinusOne << "\n"
    "\n"
    "VPD.UnderlyingPotential = VPD." << underlyingPotential << "\n"
    "\n"
    "resultString = \"error\"\n"
    "if ( ctMajorVersion is 2 ):\n"
    "    tunnelingResult = CTPD.fullTunneling( path_pts = trueAndFalseVacua,\n"
    "                                  V = VPD.PotentialForCosmotransitions,\n"
    "                                  dV = VPD.GradientForCosmotransitions,\n"
    "                                          tunneling_init_params = dict(\n"
    "                          alpha = tunnelingSymmetryDimensionMinusOne ),\n"
    "                                   tunneling_findProfile_params = dict(\n"
    "                                          npoints = tunnelPathPoints ),\n"
    "                                      deformation_deform_params = dict(\n"
    "                                  maxiter = innerLoopMaxDeformations ),\n"
    "                                   maxiter = outerLoopMaxDeformations )\n"
    "    resultString = str( tunnelingResult.action )\n"
    "elif ( ctMajorVersion is 1 ):\n"
    "    tunnelingCalculator = CTPD.fullTunneling( phi = trueAndFalseVacua,\n"
    "                                  V = VPD.PotentialForCosmotransitions,\n"
    "                                  dV = VPD.GradientForCosmotransitions,\n"
    "                            alpha = tunnelingSymmetryDimensionMinusOne,\n"
    "                                           npoints = tunnelPathPoints )\n"
    "    tunnelingCalculator.run( maxiter = innerLoopMaxDeformations,\n"
    "                             maxiter2 = outerLoopMaxDeformations )\n"
    "    resultString = str( tunnelingCalculator.findAction() )\n"
    "\n"
    "outputFile = open( \"" << pythonResultFilename << "\", \"w\" )\n"
    "outputFile.write( resultString )\n"
    "outputFile.close()\n"
    "\n"
    "# End of automatically generated file!\n";

    pythonFile.close();
    systemCommand.assign( "python " );
    systemCommand.append( pythonMainFilename );
    std::cout
    << std::endl
    << "About to run custom Python program calling CosmoTransitions!"
    << std::endl << "Unfortunately it is likely to take quite some time (at"
    << " least 10 minutes for 4 fields at 1-loop order, probably at least an"
    << " hour for 6 fields) and the output to the terminal can lag a lot (it"
    << " might only show up after the Python has finished even)."
    << std::endl << "Calling system( \" " << systemCommand << " \")..."
    << std::endl << "-----------------";
    std::cout << std::endl;
    BOL::UsefulStuff::runSystemCommand( systemCommand );
    std::cout
    << std::endl << "-----------------" << std::endl
    << "Parsing output from " << pythonMainFilename << ".";
    std::cout << std::endl;

    double calculatedAction( -1.0 );
    std::ifstream resultStream;
    resultStream.open( pythonResultFilename.c_str() );
    resultStream >> calculatedAction;
    resultStream.close();

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Action calculated by CosmoTransitions = " << calculatedAction;
    std::cout << std::endl;/**/

    return calculatedAction;
  }

} /* namespace VevaciousPlusPlus */
