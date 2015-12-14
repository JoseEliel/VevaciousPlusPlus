/*
 * VevaciousPlusPlusMain.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */


#include "CommonIncludes.hpp"
#include "VevaciousPlusPlus.hpp"
#include "LHPC/SimpleLhaParser.hpp"


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
  std::string inputFile( argumentParser.fromTag( "InputFile",
                                                "" ) );
  if( !(inputFile.empty()) )
  {
    vevaciousPlusPlus.RunPoint( inputFile );
    vevaciousPlusPlus.WriteResultsAsXmlFile( argumentParser.fromTag( "OutputFile",
                                                   ( inputFile + ".vout" ) ) );
    vevaciousPlusPlus.AppendResultsToLhaFile( inputFile );
  }

  // Solve a directory full of parameter points, if one was given with the
  // <InputFolder> tag.
  std::string inputFolder( argumentParser.fromTag( "InputFolder",
                                                   "" ) );
  std::string outputFolder( argumentParser.fromTag( "OutputFolder",
                                                    "" ) );

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
      vevaciousPlusPlus.WriteResultsAsXmlFile( placeholderManager.getOutput()
                                         + ".vout" );
      std::string copyCommand( "cp " );
      copyCommand.append( placeholderManager.getInput() );
      copyCommand.append( " " );
      copyCommand.append( placeholderManager.getOutput() );
      BOL::UsefulStuff::runSystemCommand( copyCommand );
      vevaciousPlusPlus.AppendResultsToLhaFile( placeholderManager.getOutput() );
    }
  }

  std::cout
  << std::endl
  << "Vevacious finished running.";
  std::cout << std::endl;

  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}
