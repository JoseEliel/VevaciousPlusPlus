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
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                                     double const survivalProbabilityThreshold,
                                        unsigned int const temperatureAccuracy,
                                     std::string const& pathToCosmotransitions,
                                      unsigned int const resolutionOfDsbVacuum,
                                              unsigned int const maxInnerLoops,
                                              unsigned int const maxOuterLoops,
                           unsigned int const thermalStraightPathFitResolution,
                                      double const vacuumSeparationFraction ) :
    BounceActionTunneler( tunnelingStrategy,
                          survivalProbabilityThreshold,
                          temperatureAccuracy,
                          vacuumSeparationFraction ),
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
  // integrated over three dimensions (for non-zero temperature) for
  // tunneling from falseVacuum to trueVacuum at temperature
  // tunnelingTemperature. It does so by writing and running a Python program
  // using the potential from pythonPotentialFilename for CosmoTransitions to
  // use to calculate the bounce action at tunnelingTemperature. The vacua
  // are assumed to already be the minima at tunnelingTemperature.
  double CosmoTransitionsRunner::BounceAction(
                                    PotentialFunction const& potentialFunction,
                                        PotentialMinimum const& falseVacuum,
                                        PotentialMinimum const& trueVacuum,
                                        double const tunnelingTemperature )
  {
    std::string const
    pythonResultFilename( "VevaciousCosmoTransitionsResult.txt" );
    std::string const
    pythonMainFilename( "VevaciousCosmoTransitionsRunner.py" );
    std::string systemCommand( "rm " );
    systemCommand.append( pythonMainFilename );
    systemCommand.append( "c" );
    system( systemCommand.c_str() );
    std::ofstream pythonFile;
    pythonFile.open( pythonMainFilename.c_str() );
    pythonFile << std::setprecision( 12 );
    pythonFile << "# Automatically generated file! Modify at your own PERIL!\n"
    "# This file was created by Vevacious version "
    << VersionInformation::CurrentVersion() << "\n"
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
    unsigned int tunnelingSymmetryDimensionMinusOne( 3 );
    std::string underlyingPotential( "JustLoopCorrectedPotential" );
    if( tunnelingTemperature > 0.0 )
    {
      pythonFile
      << "VPD.SetGlobalTemperature( " << tunnelingTemperature << " )\n";
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
    system( systemCommand.c_str() );
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
  // thermal dependence of the action, then uses Minuit2 to find the optimal
  // tunneling temperature, then writes and runs another Python program to
  // use CosmoTransitions to calculate the thermal action at this optimal
  // temperature.
  void CosmoTransitionsRunner::ContinueThermalTunneling(
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                              double const potentialAtOriginAtZeroTemperature )
  {
    // We are going to fit a function that diverges at
    // criticalTunnelingTemperature, so we shouldn't pick a point too close to
    // close to the divergence.
    double const highestFitTemperature( 0.9
                                * rangeOfMaxTemperatureForOriginToTrue.first );

    // Ideally the DSB vacuum evaporates at a much lower temperature than the
    // critical tunneling temperature for the true vacuum, so we can restrict
    // our fit to temperatures above the evaporation temperature where there is
    // a kink in the thermal action due to the rapid acceleration of the false
    // vacuum end of the tunneling path followed by a complete halt, with
    // respect to increasing temperature.
    std::vector::size_type
    straightPathPoints( static_cast< std::vector::size_type >(
                                          thermalStraightPathFitResolution ) );
    double
    lowestFitTemperature( 1.1 * rangeOfMaxTemperatureForOriginToFalse.second );
    if( rangeOfMaxTemperatureForOriginToTrue.first
        < ( 2.0 * rangeOfMaxTemperatureForOriginToFalse.second ) )
    {
      lowestFitTemperature = 1.0;
      if( highestFitTemperature
          > rangeOfMaxTemperatureForOriginToFalse.second )
      {
        straightPathPoints *= 2;
      }
    }
    std::vector< double > fitTemperatures( straightPathPoints );
    double const
    stepTemperature( ( highestFitTemperature - lowestFitTemperature )
                     / static_cast< double >( straightPathPoints - 1 ) );
    double currentTemperature( lowestFitTemperature );
    for( std::vector::size_type whichNode( 0 );
         whichNode < straightPathPoints;
         ++whichNode )
    {
      fitTemperatures[ whichNode ] = currentTemperature;
      currentTemperature += stepTemperature;
    }
    std::vector< double > straightPathActions;
    CalculateStraightPathActions( potentialFunction,
                                  falseVacuum,
                                  trueVacuum,
                                  fitTemperatures,
                                  straightPathActions );

    ThermalActionFitter thermalActionFitter( fitTemperatures,
                                             straightPathActions,
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
    MinuitPotentialMinimizer thermalPotentialMinimizer( potentialFunction );
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

    double thermalAction( BounceAction( potentialFunction,
                                        thermalFalseVacuum,
                  thermalPotentialMinimizer( trueVacuum.FieldConfiguration() ),
                                    dominantTemperatureInGigaElectronVolts ) );
    logOfMinusLogOfThermalProbability = ( lnOfThermalIntegrationFactor
                                          - ( thermalAction
                                     / dominantTemperatureInGigaElectronVolts )
                                          - log( thermalAction ) );
    SetThermalSurvivalProbability();
  }

  // This uses a BubbleShootingOnSpline object at different temperatures to
  // fill straightPathActions based on straight paths between the thermal
  // vacua. It uses the given potentialFunction rather than pythonPotential.
  void CosmoTransitionsRunner::InteralGuessFromStraightPaths(
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                  std::vector< double > const& fitTemperatures,
                                   std::vector< double >& straightPathActions )
  {
    // First we set up the (square of the) threshold distance that we demand
    // between the vacua at every temperature to trust the tunneling
    // calculation.
    double const thresholdSeparationSquared( vacuumSeparationFractionSquared
                                * falseVacuum.SquareDistanceTo( trueVacuum ) );

    straightPathActions.clear();
    BubbleShootingOnPathInFieldSpace
    actionCalculator( ( 1.0 / static_cast< double >( resolutionOfDsbVacuum ) ),
                      32 );
    // The number of shoot attempts doesn't need to be fixed at 32, but it is
    // unlikely that anyone will ever want to change it.
    PotentialMinimum thermalFalseVacuum( falseVacuum );
    PotentialMinimum thermalTrueVacuum( trueVacuum );
    std::vector< std::vector< double > > straightPath( 2 );
    MinuitPotentialMinimizer thermalPotentialMinimizer( potentialFunction );
    for( std::vector< double >::const_iterator
         fitTemperature( fitTemperatures.begin() );
         fitTemperature < fitTemperatures.end();
         ++fitTemperature )
    {
      thermalPotentialMinimizer.SetTemperature( *fitTemperature );
      thermalFalseVacuum
      = thermalPotentialMinimizer( thermalFalseVacuum.FieldConfiguration() );
      thermalTrueVacuum
      = thermalPotentialMinimizer( thermalTrueVacuum.FieldConfiguration() );
      if( ( thermalFalseVacuum.PotentialValue()
            <= thermalTrueVacuum.PotentialValue() )
          ||
          ( thermalTrueVacuum.SquareDistanceTo( thermalFalseVacuum )
            < thresholdSeparationSquared ) )
      {
        // If the thermal vacua have gotten so close that a tunneling
        // calculation is suspect, or tunneling to the panic vacuum has become
        // impossible, we break and take only the contributions from lower
        // temperatures.
        break;
      }
      straightPath.front() = thermalFalseVacuum.FieldConfiguration();
      straightPath.back() = thermalTrueVacuum.FieldConfiguration();
      LinearSplineThroughNodes straightSplinePath( straightPath,
                                                   std::vector< double >( 0 ),
                                                   *fitTemperature );
      SplinePotential potentialApproximation( potentialFunction,
                                              straightSplinePath,
                                              resolutionOfDsbVacuum,
                                              thresholdSeparationSquared );
      actionCalculator.ResetVacua( potentialFunction,
                                   thermalFalseVacuum,
                                   thermalTrueVacuum,
                                   *fitTemperature );
      BubbleProfile* bubbleProfile( actionCalculator( straightSplinePath,
                                                    potentialApproximation ) );
      straightPathActions.push_back( bubbleProfile->BounceAction() );
      delete bubbleProfile;
    }
  }

  // This writes a Python programme using CosmoTransitions with the minimum
  // number of deformations to get a set of actions at temperatures, and then
  // reads in the file created to fill straightPathActions.
  void CosmoTransitionsRunner::FitFromCosmoTransitionsStraightPaths(
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                  std::vector< double > const& fitTemperatures,
                                   std::vector< double >& straightPathActions )
  {
    // First we set up the (square of the) threshold distance that we demand
    // between the vacua at every temperature to trust the tunneling
    // calculation.
    double const thresholdSeparationSquared( vacuumSeparationFractionSquared
                                * falseVacuum.SquareDistanceTo( trueVacuum ) );

    std::string const
    pythonResultFilename( "VevaciousCosmoTransitionsThermalFitResult.txt" );
    std::string const
    pythonMainFilename( "VevaciousCosmoTransitionsThermalFitter.py" );
    std::string systemCommand( "rm " );
    systemCommand.append( pythonMainFilename );
    systemCommand.append( "c" );
    system( systemCommand.c_str() );
    std::ofstream pythonFile;
    pythonFile.open( pythonMainFilename.c_str() );
    pythonFile << std::setprecision( 12 );
    pythonFile << "# Automatically generated file! Modify at your own PERIL!\n"
    "# This file was created by Vevacious version "
    << VersionInformation::CurrentVersion() << "\n"
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

    MinuitPotentialMinimizer thermalPotentialMinimizer( potentialFunction );
    PotentialMinimum thermalFalseVacuum( thermalPotentialMinimizer(
                                          falseVacuum.FieldConfiguration() ) );
    PotentialMinimum thermalTrueVacuum( thermalPotentialMinimizer(
                                           trueVacuum.FieldConfiguration() ) );
    straightPathActions.clear();
    for( std::vector::size_type whichNode( 0 );
         whichNode < fitTemperatures.size();
         ++whichNode )
    {
      thermalPotentialMinimizer.SetTemperature( fitTemperatures[ whichNode ] );
      // We update the positions of the thermal vacua based on their positions
      // at the last temperature step.
      thermalFalseVacuum
      = thermalPotentialMinimizer( thermalFalseVacuum.FieldConfiguration() );
      thermalTrueVacuum
      = thermalPotentialMinimizer( thermalTrueVacuum.FieldConfiguration() );

      if( thermalTrueVacuum.SquareDistanceTo( thermalFalseVacuum )
          < thresholdSeparationSquared )
      {
        // If the thermal vacua have gotten so close that a tunneling
        // calculation is suspect, we break and take only the contributions
        // from lower temperatures.
        break;
      }
      straightPathActions.push_back( std::numeric_limits< double >::max() );
      if( whichNode > 0 )
      {
        pythonFile << ",\n ";
      }
      pythonFile << "[ " << fitTemperatures[ whichNode ] << ", [";
      std::vector< double > const&
      thermalTrueVacuumFields( thermalTrueVacuum.FieldConfiguration() );
      std::vector< double > const&
      thermalFalseVacuumFields( thermalFalseVacuum.FieldConfiguration() );
      for( std::vector< double >::const_iterator
           fieldValue( thermalTrueVacuumFields.begin() );
           fieldValue < thermalTrueVacuumFields.end();
           ++fieldValue )
      {
        if( fieldValue != thermalTrueVacuumFields.begin() )
        {
          pythonFile << ", ";
        }
        pythonFile << *fieldValue;
      }
      pythonFile << "], [";
      for( std::vector< double >::const_iterator
           fieldValue( thermalFalseVacuumFields.begin() );
           fieldValue < thermalFalseVacuumFields.end();
           ++fieldValue )
      {
        if( fieldValue != thermalFalseVacuumFields.begin() )
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
    system( systemCommand.c_str() );
    std::cout
    << std::endl << "-----------------" << std::endl
    << "Parsing output from " << pythonMainFilename << ".";
    std::cout << std::endl;

    std::ifstream resultStream;
    resultStream.open( pythonResultFilename.c_str() );
    for( std::vector::size_type whichNode( 0 );
         whichNode < straightPathActions.size();
         ++whichNode )
    {
      resultStream >> straightPathActions[ whichNode ];
    }
    resultStream.close();
  }

} /* namespace VevaciousPlusPlus */
