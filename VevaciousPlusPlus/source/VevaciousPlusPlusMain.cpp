/*
 * VevaciousPlusPlusMain.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */


#include "../include/CommonIncludes.hpp"
#include "../include/VevaciousPlusPlus.hpp"

int main( int argumentCount,
          char** argumentCharArrays )
{
  // This takes care of the command-line arguments.
  BOL::ArgumentParser argumentParser( argumentCount,
                                      argumentCharArrays,
                                      "init",
                                      "VevaciousPlusPlusInitialization.xml" );


  // To create the VevaciousPlusPlus object, one must know what kind of
  // strategies one wants to use to minimize the potential and to calculate
  // tunneling between vacua. These strategies are chosen by passing instances
  // of classes derived from abstract base classes to the VevaciousPlusPlus
  // constructor. Since constructing such instances requires one to provide the
  // code enabling them to perform their tasks, one has to start with an
  // instance of a class derived from the SlhaManager anstract base class (so
  // that SLHA data is shared properly), then pass that on to an instance of a
  // class derived from the PotentialFunction abstract base class, and then
  // pass both to the classes derived from the PotentialMinimizer and
  // TunnelingCalculator abstract base classes, unless one wants to put in some
  // hacking effort.

  // The default strategies are:
  // - minimizing the potential by a combination of homotopy continuation to
  //   find starting points with gradient minimization, as implemented in the
  //   HomotopyContinuationAndGradient class, which requires an instance of a
  //   class derived from the HomotopyContinuationReadyPotential abstract base
  //   class;
  // - calculating tunneling by evaluating the bounce action along paths
  //   between the vacua parameterized by splines, as implemented in the
  //   BounceWithSplines class, which requires an instance of a class derived
  //   from the PotentialFunction abstract base clases.
  // Hence creating an instance of the FixedScaleOneLoopPotential class,
  // which inherits from both HomotopyContinuationReadyPolynomial and
  // PotentialFunction (through HomotopyContinuationReadyPolynomial itself
  // already being derived from PotentialFunction through
  // HomotopyContinuationReadyPotential), with an instance of the
  // RunningParameterManager class, forms the basis of the components of the
  // VevaciousPlusPlus object.

  // The RunningParameterManager class handles SLHA data and interpolates to
  // the scale dependence as given in the SLHA file.
  VevaciousPlusPlus::RunningParameterManager runningParameterManager;

  // This file assumes that the default homotopy continuation followed by
  // gradient minimization followed by tunneling using the various provided
  // classes will be used. If you make custom calculators, you're probably fine
  // with modifying this file to use them.

  // The hierarchy of classes is:
  // VevaciousPlusPlus needs
  // a PotentialMinimizer and
  // a TunnelingCalculator.
  // The only provided strategy for a
  // PotentialMinimizer is the
  // HomotopyContinuationAndGradient class, which needs
  // a HomotopyContinuationSolver and
  // a GradientBasedMinimizer.
  // There are 2 choices for the HomotopyContinuationSolver:
  // 1) BasicPolynomialHomotopyContinuation: a direct implementation written by
  //                                         BOL (not yet ready!)
  // or 2) Hom4ps2Runner: a class that runs the external binary HOM4PS2 (by
  //                      Tsung-Lin Lee, Tien-Yien Li, and Chih-Hsiung Tsai).
  // Both choices for HomotopyContinuationSolver require a
  // HomotopyContinuationReadyPolynomial as the PotentialFunction.
  // There is only 1 choice for GradientBasedMinimizer:
  // 1) MinuitPotentialMinimizer, which uses the Minuit2 library by Fred James.
  // The only provided strategy for a
  // TunnelingCalculator is the
  // BounceWithSplines class.
  // There are 2 choices for the BounceWithSplines:
  // 1) MinuitBounceActionMinimizer: a class which uses Minuit2 to minimize the
  //                                 bounce action directly (not yet ready!)
  // or 2) CosmoTransitionsRunner: a class writes and runs an external Python
  //                               program to use the CosmoTransitions library
  //                               to calculate the minimal bounce action.
  // Both choices for HomotopyContinuationSolver require a PotentialFunction.
  // Hence a HomotopyContinuationReadyPolynomial must be chosen, so that the
  // HomotopyContinuationSolver and BounceWithSplines instances can be
  // constructed, which can then themselves be used to construct the
  // PotentialMinimizer and the TunnelingCalculator, and thus the
  // VevaciousPlusPlus instance.

  // The only HomotopyContinuationReadyPolynomial classes provided so far are
  // both derived from PotentialFromPolynomialAndMasses.
  VevaciousPlusPlus::PotentialFromPolynomialAndMasses*
  potentialFunction( NULL );

  std::string potentialClass( argumentParser.fromTag( "PotentialClass",
                                              "FixedScaleOneLoopPotential" ) );
  std::string modelFile( argumentParser.fromTag( "model",
                                                 "./ModelFiles/SM.vin" ) );

  if( potentialClass.compare( "FixedScaleOneLoopPotential" ) == 0 )
  {
    potentialFunction
    = new VevaciousPlusPlus::FixedScaleOneLoopPotential( modelFile,
                                                     runningParameterManager );

    std::cout
    << std::endl
    << "Created FixedScaleOneLoopPotential from " << modelFile;
    std::cout << std::endl;
  }
  else if( potentialClass.compare( "RgeImprovedOneLoopPotential" ) == 0 )
  {
    potentialFunction
    = new VevaciousPlusPlus::RgeImprovedOneLoopPotential( modelFile,
                                                     runningParameterManager );

    std::cout
    << std::endl
    << "Created RgeImprovedOneLoopPotential from " << modelFile;
    std::cout << std::endl;
  }
  else
  {
    std::cout
    << std::endl
    << "PotentialType was not a recognized form! The only currently-valid"
    << " types are \"FixedScaleOneLoopPotential\" and"
    << " \"RgeImprovedOneLoopPotential\". Aborting!";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }


  VevaciousPlusPlus::HomotopyContinuationSolver*
  homotopyContinuationSolver( NULL );

  std::string homotopyContinuationClass(
                           argumentParser.fromTag( "HomotopyContinuationClass",
                                     "BasicPolynomialHomotopyContinuation" ) );
  if( homotopyContinuationClass.compare( "FixedScaleOneLoopPotential" ) == 0 )
  {
    homotopyContinuationSolver
    = new VevaciousPlusPlus::BasicPolynomialHomotopyContinuation(
                       potentialFunction->HomotopyContinuationTargetSystem() );

    std::cout
    << std::endl
    << "Created BasicPolynomialHomotopyContinuation, but this does not yet"
    << " work...";
    std::cout << std::endl;
  }
  else if( homotopyContinuationClass.compare( "Hom4ps2Runner" ) == 0 )
  {
    std::string pathToHom4ps2( argumentParser.fromTag( "PathToHom4ps2",
                                                              "./HOM4PS2/" ) );
    homotopyContinuationSolver
    = new VevaciousPlusPlus::Hom4ps2Runner(
                         potentialFunction->HomotopyContinuationTargetSystem(),
                                            pathToHom4ps2,
                                     argumentParser.fromTag( "Hom4ps2Argument",
                                                               "1" ) );

    std::cout
    << std::endl
    << "Created Hom4ps2Runner to run " << pathToHom4ps2 << "/hom4ps2";
    std::cout << std::endl;
  }
  else
  {
    std::cout
    << std::endl
    << "HomotopyContinuationClass was not a recognized form! The only"
    << " currently-valid types are \"BasicPolynomialHomotopyContinuation\""
    << " (not yet implemented, sorry!) and \"Hom4ps2Runner\". Aborting!";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }

  // Now the HomotopyContinuationAndGradient object can be constructed:
  VevaciousPlusPlus::HomotopyContinuationAndMinuit
  potentialMinimizer( *potentialFunction,
                      *homotopyContinuationSolver,
                     BOL::StringParser::stringToDouble( argumentParser.fromTag(
                                                   "MinimaSeparationThreshold",
                                                                   "0.1" ) ) );


  std::string
  tunnelingStrategyAsString( argumentParser.fromTag( "TunnelingStrategy",
                                                     "DefaultTunneling" ) );
  VevaciousPlusPlus::TunnelingCalculator::TunnelingStrategy
  tunnelingStrategy( VevaciousPlusPlus::TunnelingCalculator::NotSet );
  if( ( tunnelingStrategyAsString.compare( "DefaultTunneling" ) == 0 )
      ||
      ( tunnelingStrategyAsString.compare( "ThermalThenQuantum" ) == 0 ) )
  {
    tunnelingStrategy
    = VevaciousPlusPlus::TunnelingCalculator::ThermalThenQuantum;
  }
  else if( tunnelingStrategyAsString.compare( "QuantumThenThermal" ) == 0 )
  {
    tunnelingStrategy
    = VevaciousPlusPlus::TunnelingCalculator::QuantumThenThermal;
  }
  else if( tunnelingStrategyAsString.compare( "JustThermal" ) == 0 )
  {
    tunnelingStrategy
    = VevaciousPlusPlus::TunnelingCalculator::JustThermal;
  }
  else if( tunnelingStrategyAsString.compare( "JustQuantum" ) == 0 )
  {
    tunnelingStrategy
    = VevaciousPlusPlus::TunnelingCalculator::JustQuantum;
  }
  else if( ( tunnelingStrategyAsString.compare( "NoTunneling" ) == 0 )
           ||
           ( tunnelingStrategyAsString.compare( "None" ) == 0 ) )
  {
    tunnelingStrategy
    = VevaciousPlusPlus::TunnelingCalculator::NoTunneling;
  }
  else
  {
    std::cout
    << std::endl
    << "TunnelingStrategy was not a recognized form! The only currently-valid"
    << " types are \"DefaultTunneling\" (=\"ThermalThenQuantum\"),"
    << " \"ThermalThenQuantum\", \"QuantumThenThermal\", \"JustThermal\","
    << " \"JustQuantum\", \"NoTunneling\", and \"None\" ( =\"NoTunneling\")."
    << " Aborting!";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }
  std::cout
  << std::endl
  << "Tunneling strategy is " << tunnelingStrategyAsString;
  std::cout << std::endl;

  std::string survivalProbabilityThresholdAsString(
                        argumentParser.fromTag( "SurvivalProbabilityThreshold",
                                                "0.01" ) );
  double survivalProbabilityThreshold( -1.0 );
  bool validInput( BOL::StringParser::stringIsDouble(
                                          survivalProbabilityThresholdAsString,
                                              survivalProbabilityThreshold ) );
  if( !validInput
      ||
      ( survivalProbabilityThreshold <= 0.0 )
      ||
      ( survivalProbabilityThreshold >= 1.0 ) )
  {
    std::cout
    << std::endl
    << "SurvivalProbabilityThreshold was not a number between 0.0 and 1.0 but"
    << " was " << survivalProbabilityThresholdAsString << "."
    << " Aborting!";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }


  VevaciousPlusPlus::TunnelingCalculator*
  tunnelingCalculator( NULL );

  std::string tunnelingClass( argumentParser.fromTag( "TunnelingClass",
                                             "MinuitBounceActionMinimizer" ) );

  if( tunnelingClass.compare( "MinuitBounceActionMinimizer" ) == 0 )
  {
    tunnelingCalculator
    = new VevaciousPlusPlus::MinuitBounceActionMinimizer( *potentialFunction,
                                                          tunnelingStrategy,
                                                survivalProbabilityThreshold );

    std::cout
    << std::endl
    << "Created MinuitBounceActionMinimizer, but this does not yet"
    << " work...";
    std::cout << std::endl;
    std::cout << std::endl;
  }
  else if( tunnelingClass.compare( "CosmoTransitionsRunner" ) == 0 )
  {
    std::string
    pathToCosmotransitions( argumentParser.fromTag( "PathToCosmotransitions",
                                               "./CosmoTransitions-1.0.2/" ) );
    tunnelingCalculator
    = new VevaciousPlusPlus::CosmoTransitionsRunner( *potentialFunction,
                                                     *potentialFunction,
                                                     tunnelingStrategy,
                                                  survivalProbabilityThreshold,
                                                     pathToCosmotransitions );

    std::cout
    << std::endl
    << "Created CosmoTransitionsRunner from " << modelFile;
    std::cout << std::endl;
  }
  else
  {
    std::cout
    << std::endl
    << "TunnelingClass was not a recognized form! The only currently-valid"
    << " types are \"MinuitBounceActionMinimizer\" (not yet implemented,"
    << " sorry!)  and \"CosmoTransitionsRunner\". Aborting!";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }


  // Create the VevaciousTwo object, telling it where to find its settings,
  // such as the model to use and the way to calculate the tunneling time.
  VevaciousPlusPlus::VevaciousPlusPlus vevaciousPlusPlus( argumentParser,
                                                       runningParameterManager,
                                                          potentialMinimizer,
                                                        *tunnelingCalculator );

  // Solve a parameter point, if one was given directly with the <slha> tag:
  std::string slhaFile( argumentParser.fromTag( "slha",
                                                "" ) );
  if( !(slhaFile.empty()) )
  {
    vevaciousPlusPlus.RunPoint( slhaFile );
    vevaciousPlusPlus.WriteXmlResults( argumentParser.fromTag( "output",
                                                    ( slhaFile + ".vout" ) ) );
    vevaciousPlusPlus.WriteSlhaResults( slhaFile );
  }

  // Solve a directory full of parameter points, if one was given with the
  // <InputFolder> tag.
  std::string inputFolder(  argumentParser.fromTag( "InputFolder",
                                                       "" ) );
  std::string outputFolder(  argumentParser.fromTag( "OutputFolder",
                                                        "" ) );
  // debugging:
  /**/std::cout << std::endl << "debugging:"
  << std::endl
  << "inputFolder = \"" << inputFolder << "\", inputFolder.empty() = "
  << inputFolder.empty();
  std::cout << std::endl;/**/

  if( !(inputFolder.empty()) )
  {
    if( outputFolder.empty() )
    {
      std::cout
      << std::endl
      << "OutputFolder string must not be empty string! Use \"./\" for the"
      << " current working folder.";
      std::cout << std::endl;

      return EXIT_FAILURE;
    }
    if( outputFolder.compare( inputFolder ) == 0 )
    {
      std::cout
      << std::endl
      << "Input folder and output folder must be different, as the filenames"
      << " do not change!";
      std::cout << std::endl;

      return EXIT_FAILURE;
    }
    BOL::FilePlaceholderManager placeholderManager( "",
                                                    ".placeholder",
                                                    "" );
    placeholderManager.prepareFilenames( inputFolder,
                                         outputFolder,
                                         outputFolder );

    while( placeholderManager.holdNextPlace() )
    {
      vevaciousPlusPlus.RunPoint( placeholderManager.getInput() );
      vevaciousPlusPlus.WriteXmlResults( placeholderManager.getOutput()
                                         + ".vout" );
      std::string copyCommand( "cp " );
      copyCommand.append( placeholderManager.getInput() );
      copyCommand.append( " " );
      copyCommand.append( placeholderManager.getOutput() );
      BOL::UsefulStuff::runSystemCommand( copyCommand );
      vevaciousPlusPlus.WriteSlhaResults( placeholderManager.getOutput() );
    }
  }


  // The FixedScaleOneLoopPotential constructor takes a string with the name of
  // the model file (including the path) and a reference to a
  // RunningParameterManager instance:
  VevaciousPlusPlus::FixedScaleOneLoopPotential
  fixedScalePotential( *potentialFunction );

  // Both the FixedScaleOneLoopPotential and the RgeImprovedOneLoopPotential
  // classes can take a reference to a PotentialFromPolynomialAndMasses
  // instance to make a copy constructor, and since they are both derived from
  // PotentialFromPolynomialAndMasses, we can make an instance of 1 of them
  // from an instance of the other.
  VevaciousPlusPlus::RgeImprovedOneLoopPotential
  rgeImprovedPotential( *potentialFunction );

  // debugging:
  /**/std::cout << std::endl << "debugging:"
  << std::endl;
  runningParameterManager.UpdateSlhaData( slhaFile );

  std::cout
  << "Fixed-scale at origin = "
  << fixedScalePotential( fixedScalePotential.FieldValuesOrigin() )
  << std::endl;
  std::cout
  << "Fixed-scale at DSB = "
  << fixedScalePotential( fixedScalePotential.DsbFieldValues() )
  << std::endl;
  /**/

  std::cout
  << std::endl
  << "------" << std::endl << std::endl;
  std::cout << std::endl;

  double subtractionConstant( fixedScalePotential(
                                   fixedScalePotential.FieldValuesOrigin() ) );
  std::cout
  << std::endl
  << "At field origin, fixedScalePotential = " << subtractionConstant;
  std::cout << std::endl;


  std::vector< double >
  testConfiguration( fixedScalePotential.FieldValuesOrigin() );
  for( int vdStep( 0 );
       vdStep < 13;
       ++vdStep )
  {
    testConfiguration[ 0 ] = ( 0.1
                               * (double)vdStep
                               //* fixedScalePotential.DsbFieldValues()[ 0 ] );
                               * 1.0 );
    for( int vuStep( 0 );
         vuStep < 13;
         ++vuStep )
    {
      testConfiguration[ 1 ] = ( 0.1
                                 * (double)vuStep
                                 //* fixedScalePotential.DsbFieldValues()[ 1 ] );
                                 * 1.0 );
      std::cout << "For "
      << fixedScalePotential.FieldConfigurationAsMathematica(
                                                            testConfiguration )
      << ", fixedScalePotential => "
      << fixedScalePotential( testConfiguration )
      << " => " << ( fixedScalePotential( testConfiguration )
                     - subtractionConstant )
      << "; tree = "
      << fixedScalePotential.QuickApproximation( testConfiguration );
      std::cout << std::endl;
    }
  }/**/

  std::cout << std::endl;
  std::cout << std::endl;

  delete tunnelingCalculator;
  delete homotopyContinuationSolver;
  delete potentialFunction;

  std::cout
  << std::endl
  << "Still to do:" << std::endl
  << "write MinuitBounceActionMinimizer - remember to check action after e.g."
  << " 10 * n_fields FCN calls" << std::endl
  << "add in <TakenPositive>" << std::endl
  << "add in <RollOnlyTreeLevelMinima>" << std::endl
  << "write BasicPolynomialHomotopyContinuation" << std::endl
  << "think about uncertainties" << std::endl;
  std::cout << std::endl;


  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}






