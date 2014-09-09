/*
 * VevaciousPlusPlusMain.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */


#include "CommonIncludes.hpp"
#include "VevaciousPlusPlus.hpp"

int main( int argumentCount,
          char** argumentCharArrays )
{
  // This takes care of the command-line arguments.
  BOL::ArgumentParser argumentParser( argumentCount,
                                      argumentCharArrays,
                                      "input",
                                      "VevaciousPlusPlusMainInput.xml" );

  // The default is to construct the VevaciousPlusPlus object with an input
  // initialization file name, which will create a PotentialMinimizer and
  // TunnelingCalculator internal to the VevaciousPlusPlus object, managing
  // the memory allocated in doing so, with all components initialized
  // according to the XML input file. Alternatively, one can create the
  // PotentialMinimizer and TunnelingCalculator components externally and pass
  // them to the other constructor.
  VevaciousPlusPlus::VevaciousPlusPlus
  vevaciousPlusPlus( argumentParser.fromTag( "InitializationFile",
  "./InitializationFiles/VevaciousPlusPlusDefaultObjectInitialization.xml" ) );

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
  std::string inputFolder( argumentParser.fromTag( "InputFolder",
                                                   "" ) );
  std::string outputFolder( argumentParser.fromTag( "OutputFolder",
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


  // placeholder:
  /**/std::cout << std::endl
  << "Placeholder: "
  << "EVERYTHING AFTER THIS IS CODE TO TEST STUFF!";
  std::cout << std::endl;/**/

  // The RunningParameterManager class handles SLHA data and interpolates to
  // the scale dependence as given in the SLHA file.
  VevaciousPlusPlus::RunningParameterManager runningParameterManager;


  VevaciousPlusPlus::FixedScaleOneLoopPotential
  fixedScalePotential( argumentParser.fromTag( "model",
                               "/Users/oleary/BOL/Cplusplus/VevaciousPlusPlus/"
                     "VevaciousPlusPlus/ModelFiles/RealMssmWithStopVevs.vin" ),
                       10.0,
                       true,
                       runningParameterManager );
  VevaciousPlusPlus::RgeImprovedOneLoopPotential
  rgeImprovedPotential( fixedScalePotential );

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
    testConfiguration[ 0 ]
    = ( 0.1 * vdStep * std::max( 1.0,
                                 fixedScalePotential.DsbFieldValues()[ 0 ] ) );
    if( testConfiguration.size() > 1 )
    {
      for( int vuStep( 0 );
           vuStep < 13;
           ++vuStep )
      {
        testConfiguration[ 1 ]
        = ( 0.1 * vuStep * std::max( 1.0,
                                 fixedScalePotential.DsbFieldValues()[ 1 ] ) );
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
    }
    else
    {
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

  std::cout
  << std::endl
  << "Still to do:" << std::endl
  << "abandon path from minimizing potential on hemispheres?" << std::endl
  << "write path through nodes on bisections" << std::endl
  << "write BasicPolynomialHomotopyContinuation" << std::endl
  << "think about uncertainties" << std::endl;
  std::cout << std::endl;


  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}






