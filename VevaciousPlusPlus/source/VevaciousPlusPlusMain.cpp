/*
 * VevaciousPlusPlusMain.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */


#include "../include/StandardIncludes.hpp"
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

  VevaciousPlusPlus::HomotopyContinuationReadyPolynomial*
  potentialFunction( NULL );

  std::string potentialType( argumentParser.fromTag( "PotentialType",
                                              "FixedScaleOneLoopPotential" ) );
  if( potentialType.compare( "FixedScaleOneLoopPotential" ) == 0 )
  {
    potentialFunction = new VevaciousPlusPlus::FixedScaleOneLoopPotential(
                                               argumentParser.fromTag( "model",
                                                       "./ModelFiles/SM.vin" ),
                                                     runningParameterManager );
  }
  else if( potentialType.compare( "RgeImprovedOneLoopPotential" ) == 0 )
  {
    potentialFunction = new VevaciousPlusPlus::RgeImprovedOneLoopPotential(
                                               argumentParser.fromTag( "model",
                                                       "./ModelFiles/SM.vin" ),
                                                     runningParameterManager );
  }
  else
  {
    std::cout
    << std::endl
    << "PotentialType was not a recognized form! The only currently-valid"
    << " types are \"FixedScaleOneLoopPotential\" and"
    << " \"RgeImprovedOneLoopPotential\".";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }


  VevaciousPlusPlus::HomotopyContinuationSolver*
  homotopyContinuationSolver( NULL );

  std::string
  homotopyContinuationType( argumentParser.fromTag( "HomotopyContinuation",
                                     "BasicPolynomialHomotopyContinuation" ) );
  if( homotopyContinuationType.compare( "FixedScaleOneLoopPotential" ) == 0 )
  {
    homotopyContinuationSolver
    = new VevaciousPlusPlus::BasicPolynomialHomotopyContinuation(
                                                          *potentialFunction );
  }
  else if( homotopyContinuationType.compare( "Hom4ps2Runner" ) == 0 )
  {
    homotopyContinuationSolver
    = new VevaciousPlusPlus::Hom4ps2Runner( *potentialFunction,
                                       argumentParser.fromTag( "PathToHom4ps2",
                                                              "./HOM4PS2/" ) );
  }
  else
  {
    std::cout
    << std::endl
    << "HomotopyContinuation was not a recognized form! The only"
    << " currently-valid types are \"BasicPolynomialHomotopyContinuation\" and"
    << " \"Hom4ps2Runner\".";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }

  // Now the HomotopyContinuationAndGradient object can be constructed:
  VevaciousPlusPlus::HomotopyContinuationAndMinuit
  potentialMinimizer( *potentialFunction,
                      *homotopyContinuationSolver,
                      0.1 );


  // Also the BounceWithSplines object can now be constructed:
  VevaciousPlusPlus::BounceWithSplines
  tunnelingCalculator( *potentialFunction );

  // Create the VevaciousTwo object, telling it where to find its settings,
  // such as the model to use and the way to calculate the tunneling time.
  VevaciousPlusPlus::VevaciousPlusPlus vevaciousPlusPlus( argumentParser,
                                                       runningParameterManager,
                                                          potentialMinimizer,
                                                         tunnelingCalculator );

  // Solve a parameter point:
  std::string slhaFile( argumentParser.fromTag( "slha",
                                                "slha.out" ) );
  vevaciousPlusPlus.RunPoint( slhaFile );

  // debugging:
  /**/std::cout << std::endl << "debugging:"
  << std::endl
  << "Ran vevaciousPlusPlus.RunPoint( \"" << slhaFile << "\" )";
  std::cout << std::endl;/**/



  // The FixedScaleOneLoopPotential constructor takes a string with the name of
  // the model file (including the path) and a reference to a
  // RunningParameterManager instance:
  VevaciousPlusPlus::FixedScaleOneLoopPotential
  fixedScalePotential( argumentParser.fromTag( "model",
                                               "./ModelFiles/SM.vin" ),
                       runningParameterManager );

  // Both the FixedScaleOneLoopPotential and the RgeImprovedOneLoopPotential
  // classes can take a reference to a PotentialFromPolynomialAndMasses
  // instance to make a copy constructor, and since they are both derived from
  // PotentialFromPolynomialAndMasses, we can make an instance of 1 of them
  // from an instance of the other.
  VevaciousPlusPlus::RgeImprovedOneLoopPotential
  rgeImprovedPotential( fixedScalePotential );

  // debugging:
  /**/std::cout << std::endl << "debugging:"
  << std::endl;
  std::vector< double >
  testConfiguration( rgeImprovedPotential.NumberOfFieldVariables(),
                     0.0 );
  if( testConfiguration.size() < 2 )
  {
    testConfiguration.resize( 2,
                              0.0 );
  }
  runningParameterManager.UpdateSlhaData( slhaFile );
  rgeImprovedPotential.UpdateParameters();
  fixedScalePotential.UpdateParameters();

  std::cout
  << "Fixed-scale at origin = "
  << fixedScalePotential( fixedScalePotential.FieldValuesOrigin() )
  << std::endl;
  std::cout
  << "Fixed-scale at DSB = "
  << fixedScalePotential( fixedScalePotential.DsbFieldValues() )
  << std::endl;

  std::cout
  << std::endl
  << "------" << std::endl << std::endl;
  std::cout << std::endl;


  std::cout << "For { ";
  for( unsigned int fieldIndex( 0 );
       fieldIndex < rgeImprovedPotential.NumberOfFieldVariables();
       ++fieldIndex )
  {
    if( fieldIndex > 0 )
    {
      std::cout << ", ";
    }
    std::cout << rgeImprovedPotential.FieldName( fieldIndex ) << " -> "
    << testConfiguration[ fieldIndex ];
  }
  std::cout
  << " }, rgeImprovedPotential => "
  << rgeImprovedPotential( testConfiguration )
  << ", fixedScalePotential => "
  << fixedScalePotential( testConfiguration );
  std::cout << std::endl;
  double testTemperature( 10.0 );
  std::cout
  << "at temperature " << testTemperature << ", rgeImprovedPotential = "
  << rgeImprovedPotential( testConfiguration,
                           testTemperature ) << ", fixedScalePotential = "
  << fixedScalePotential( testConfiguration,
                          testTemperature );
  std::cout << std::endl;

  double subtractionConstant( fixedScalePotential(
                                   fixedScalePotential.FieldValuesOrigin() ) );
  std::cout
  << std::endl
  << "At field origin, fixedScalePotential = " << subtractionConstant;
  std::cout << std::endl;


  for( int vdStep( 0 );
       vdStep < 13;
       ++vdStep )
  {
    testConfiguration[ 0 ] = ( 0.1
                               * (double)vdStep
                               * fixedScalePotential.DsbFieldValues()[ 0 ] );
    for( int vuStep( 0 );
         vuStep < 13;
         ++vuStep )
    {
      testConfiguration[ 1 ] = ( 0.1
                                 * (double)vuStep
                                 * fixedScalePotential.DsbFieldValues()[ 1 ] );
      std::cout << "For { ";
      for( unsigned int fieldIndex( 0 );
           fieldIndex < fixedScalePotential.NumberOfFieldVariables();
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout
        << fixedScalePotential.FieldName( fieldIndex ) << " -> "
        << testConfiguration[ fieldIndex ];
      }
      std::cout
      << ", fixedScalePotential => "
      << fixedScalePotential( testConfiguration )
      << " => " << ( fixedScalePotential( testConfiguration )
                     - subtractionConstant )
      << "; tree = "
      << fixedScalePotential.QuickApproximation( testConfiguration );
      std::cout << std::endl;
    }
  }

  testConfiguration[ 0 ] = 300.0;
  runningParameterManager.UpdateSlhaData( slhaFile );
  rgeImprovedPotential.UpdateParameters();
  fixedScalePotential.UpdateParameters();
  std::cout << "For { ";
  for( unsigned int fieldIndex( 0 );
       fieldIndex < rgeImprovedPotential.NumberOfFieldVariables();
       ++fieldIndex )
  {
    if( fieldIndex > 0 )
    {
      std::cout << ", ";
    }
    std::cout << rgeImprovedPotential.FieldName( fieldIndex ) << " -> "
    << testConfiguration[ fieldIndex ];
  }
  std::cout
  << " }, rgeImprovedPotential => "
  << rgeImprovedPotential( testConfiguration )
  << ", fixedScalePotential => "
  << fixedScalePotential( testConfiguration );
  std::cout << std::endl;
  testTemperature = 145.0;
  std::cout
  << "at temperature " << testTemperature << ", rgeImprovedPotential = "
  << rgeImprovedPotential( testConfiguration,
                           testTemperature ) << ", fixedScalePotential = "
  << fixedScalePotential( testConfiguration,
                          testTemperature );
  std::cout << std::endl;
  runningParameterManager.UpdateSlhaData( slhaFile );
  rgeImprovedPotential.UpdateParameters();
  fixedScalePotential.UpdateParameters();
  std::cout << "For { ";
  for( unsigned int fieldIndex( 0 );
       fieldIndex < rgeImprovedPotential.NumberOfFieldVariables();
       ++fieldIndex )
  {
    if( fieldIndex > 0 )
    {
      std::cout << ", ";
    }
    std::cout << rgeImprovedPotential.FieldName( fieldIndex ) << " -> "
    << testConfiguration[ fieldIndex ];
  }
  std::cout
  << " }, rgeImprovedPotential => "
  << rgeImprovedPotential( testConfiguration )
  << ", fixedScalePotential => "
  << fixedScalePotential( testConfiguration );
  std::cout << std::endl;
  testTemperature = 205.0;
  std::cout
  << "at temperature " << testTemperature << ", rgeImprovedPotential = "
  << rgeImprovedPotential( testConfiguration,
                           testTemperature ) << ", fixedScalePotential = "
  << fixedScalePotential( testConfiguration,
                          testTemperature );
  std::cout << std::endl;
  std::cout << std::endl;

  // Write the results:
  vevaciousPlusPlus.WriteXmlResults( argumentParser.fromTag( "output",
                                                    ( slhaFile + ".vout" ) ) );
  vevaciousPlusPlus.WriteSlhaResults( slhaFile );

  delete potentialFunction;
  delete homotopyContinuationSolver;
  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}






