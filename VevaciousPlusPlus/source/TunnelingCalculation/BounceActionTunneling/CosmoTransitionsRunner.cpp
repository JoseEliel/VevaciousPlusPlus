/*
 * CosmoTransitionsRunner.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/CosmoTransitionsRunner.hpp"

namespace VevaciousPlusPlus
{
  std::string
  CosmoTransitionsRunner::pythonPotentialFilenameBase( "VevaciousPotential" );

  CosmoTransitionsRunner::CosmoTransitionsRunner(
                                       IWritesPythonPotential& pythonPotential,
                                          PotentialFunction& potentialFunction,
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                                     double const survivalProbabilityThreshold,
                                              size_t const temperatureAccuracy,
                                     std::string const& pathToCosmotransitions,
                                            size_t const resolutionOfDsbVacuum,
                                                    size_t const maxInnerLoops,
                                                    size_t const maxOuterLoops,
                              size_t const thermalStraightPathFitResolution ) :
    BounceActionTunneler( potentialFunction,
                          tunnelingStrategy,
                          survivalProbabilityThreshold,
                          temperatureAccuracy ),
    pythonPotential( pythonPotential ),
    pathToCosmotransitions( pathToCosmotransitions ),
    resolutionOfDsbVacuum( resolutionOfDsbVacuum ),
    maxInnerLoops( maxInnerLoops ),
    maxOuterLoops( maxOuterLoops ),
    thermalStraightPathFitResolution( thermalStraightPathFitResolution )
  {
    // This constructor is just an initialization list.
  }

  CosmoTransitionsRunner::~CosmoTransitionsRunner()
  {
    // This does nothing.
  }


  // This returns either the dimensionless bounce action integrated over four
  // dimensions (for zero temperature) or the dimensionful bounce action
  // integrated over three dimensions (for non-zero temperature) for tunneling
  // from falseVacuum to trueVacuum at temperature tunnelingTemperature. It
  // does so by writing and running a Python program using the potential from
  // pythonPotentialFilename for CosmoTransitions to use to calculate the
  // bounce action at tunnelingTemperature. The vacua are assumed to already be
  // the minima at tunnelingTemperature.
  double
  CosmoTransitionsRunner::BounceAction( PotentialMinimum const& falseVacuum,
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
    << VersionInformation::currentVersion << "\n"
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
    "tunnelPathPoints = " << resolutionOfDsbVacuum << "\n"
    "\n"
    "# CosmoTransitions nests a loop within a loop to do its deformations:\n"
    "innerLoopMaxDeformations = " << maxInnerLoops << "\n"
    "outerLoopMaxDeformations = " << maxOuterLoops << "\n";
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
    << std::endl << "Calling system( \"" << systemCommand << "\" )..."
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

    std::cout << std::endl << "CosmoTransitions calculated an action of "
    << calculatedAction;
    if( tunnelingTemperature > 0.0 )
    {
      std::cout << " GeV";
    }
    std::cout << "." << std::endl;

    return calculatedAction;
  }

  // This calculates the evaporation and critical temperatures, then writes
  // and runs a Python program using the potential from
  // pythonPotentialFilename for CosmoTransitions to get an estimate of the
  // thermal dependence of the action, then uses Minuit2 to find the
  // optimal tunneling temperature, then writes and runs another Python
  // program to use CosmoTransitions to calculate the thermal action at this
  // optimal temperature.
  void CosmoTransitionsRunner::ContinueThermalTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum,
                              double const potentialAtOriginAtZeroTemperature )
  {
    // We are going to fit a function that diverges at
    // criticalTunnelingTemperature, so we shouldn't pick a point too close to
    // close to the divergence.
    double const highestFitTemperature( 0.9
                                * rangeOfMaxTemperatureForOriginToTrue.first );

    std::vector< std::vector< double > >
    fitFalseVacua( thermalStraightPathFitResolution );
    std::vector< std::vector< double > >
    fitTrueVacua( thermalStraightPathFitResolution );
    std::vector< double > fitTemperatures( thermalStraightPathFitResolution );

    // Ideally the DSB vacuum evaporates at a much lower temperature than the
    // critical tunneling temperature for the true vacuum, so we can restrict
    // our fit to temperatures above the evaporation temperature where there is
    // a kink in the thermal action due to the rapid acceleration of the false
    // vacuum end of the tunneling path followed by a complete halt, with
    // respect to increasing temperature.
    if( rangeOfMaxTemperatureForOriginToTrue.first
        > ( 2.0 * rangeOfMaxTemperatureForOriginToFalse.second ) )
    {
      double const lowestFitTemperature( 1.1
                              * rangeOfMaxTemperatureForOriginToFalse.second );
      double const
      stepTemperature( ( highestFitTemperature - lowestFitTemperature )
             / static_cast< double >( thermalStraightPathFitResolution - 1 ) );
      double currentTemperature( lowestFitTemperature );
      for( size_t whichNode( 0 );
           whichNode < thermalStraightPathFitResolution;
           ++whichNode )
      {
        fitTemperatures[ whichNode ] = currentTemperature;
        thermalPotentialMinimizer.SetTemperature( currentTemperature );
        fitFalseVacua[ whichNode ] = potentialFunction.FieldValuesOrigin();
        fitTrueVacua[ whichNode ] = thermalPotentialMinimizer(
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
             / static_cast< double >( thermalStraightPathFitResolution - 1 ) );
      double currentTemperature( lowestFitTemperature );
      for( size_t whichNode( 0 );
           whichNode < thermalStraightPathFitResolution;
           ++whichNode )
      {
        fitTemperatures[ whichNode ] = currentTemperature;
        thermalPotentialMinimizer.SetTemperature( currentTemperature );
        if( currentTemperature > rangeOfMaxTemperatureForOriginToFalse.second )
        {
          fitFalseVacua[ whichNode ] = potentialFunction.FieldValuesOrigin();
        }
        else
        {
          fitFalseVacua[ whichNode ] = thermalPotentialMinimizer(
                           falseVacuum.FieldConfiguration() ).VariableValues();
        }
        fitTrueVacua[ whichNode ] = thermalPotentialMinimizer(
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
    << VersionInformation::currentVersion << "\n"
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
    "tunnelPathPoints = " << resolutionOfDsbVacuum << "\n"
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
    for( size_t whichNode( 0 );
         whichNode < thermalStraightPathFitResolution;
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
    << std::endl << "Calling system( \"" << systemCommand << "\" )..."
    << std::endl << "-----------------";
    std::cout << std::endl;
    BOL::UsefulStuff::runSystemCommand( systemCommand );
    std::cout
    << std::endl << "-----------------" << std::endl
    << "Parsing output from " << pythonMainFilename << ".";
    std::cout << std::endl;

    std::vector< double > fittedActions( thermalStraightPathFitResolution );
    std::ifstream resultStream;
    resultStream.open( pythonResultFilename.c_str() );
    for( size_t whichNode( 0 );
         whichNode < thermalStraightPathFitResolution;
         ++whichNode )
    {
      resultStream >> fittedActions[ whichNode ];
    }
    resultStream.close();

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Fitted actions as [ T, {DSB}, {panic}, S ]:" << std::endl;
    for( size_t whichNode( 0 );
         whichNode < thermalStraightPathFitResolution;
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
      std::cout << " }, { ";
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
      std::cout << " }, " << fittedActions[ whichNode ] << " ]" << std::endl;
      std::cout << std::endl;
    }
    std::cout << std::endl;/**/

    ThermalActionFitter thermalActionFitter( fitTemperatures,
                                             fittedActions,
                                  rangeOfMaxTemperatureForOriginToTrue.first );
    ROOT::Minuit2::MnMigrad
    fittedThermalActionMinimizer( thermalActionFitter,
                                  std::vector< double >( 1,
                        ( 0.5 * rangeOfMaxTemperatureForOriginToTrue.first ) ),
                                  std::vector< double >( 1,
                       ( 0.05 * rangeOfMaxTemperatureForOriginToTrue.first ) ),
                                  1 );
    dominantTemperatureInGigaElectronVolts
    = fittedThermalActionMinimizer().UserParameters().Value( 0 );

    std::cout
    << std::endl
    << "Dominant temperature for tunneling estimated to be "
    << dominantTemperatureInGigaElectronVolts << " GeV.";
    std::cout << std::endl;

    // Finally we allow CosmoTransitions to calculate the action at our best
    // guess of the optimal tunneling temperature with full path deformation.
    thermalPotentialMinimizer.SetTemperature(
                                      dominantTemperatureInGigaElectronVolts );
    // We assume that dominantTemperatureInGigaElectronVolts is high enough
    // that tunneling will be out of the field origin, as the DSB vacuum proper
    // has evaporated at this temperature.
    PotentialMinimum thermalFalseVacuum( potentialFunction.FieldValuesOrigin(),
                      potentialFunction( potentialFunction.FieldValuesOrigin(),
                                       dominantTemperatureInGigaElectronVolts )
                                - thermalPotentialMinimizer.FunctionOffset() );
    if( dominantTemperatureInGigaElectronVolts
        < rangeOfMaxTemperatureForOriginToFalse.second )
    {
      // If the DSB vacuum has not evaporated, we find out where it is exactly
      // so that it can be tunneled out of.
      thermalFalseVacuum
      = thermalPotentialMinimizer( falseVacuum.FieldConfiguration() );
    }

    double thermalAction( BounceAction( thermalFalseVacuum,
                  thermalPotentialMinimizer( trueVacuum.FieldConfiguration() ),
                                    dominantTemperatureInGigaElectronVolts ) );
    logOfMinusLogOfThermalProbability = ( lnOfThermalIntegrationFactor
                                          - ( thermalAction
                                     / dominantTemperatureInGigaElectronVolts )
                                          - log( thermalAction ) );
    SetThermalSurvivalProbability();
  }

} /* namespace VevaciousPlusPlus */
