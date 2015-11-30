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

  // debugging:
  /**/std::cout << std::endl << "debugging:"
  << std::endl
  << "Testing old and new LHA parsers.";
  std::cout << std::endl;/**/

  // std::string const slhaFileName( "CMSSM_CCB.slha.out" );
  std::string const slhaFileName( "test.slha" );

  BOL::AsciiXmlParser lhaManagerParser;
  bool const successfulOpening( lhaManagerParser.openRootElementOfFile(
      "/home/bol/BOL/C++Projects/VevaciousPlusPlus/VevaciousPlusPlus/"
      "ModelFiles/LagrangianParameters/"
      "MssmCompatibleWithSlhaOneAndSlhaTwoAndSarahOutputs.xml" ) );
  std::cout
  << std::endl
  << "successfulOpening = " << successfulOpening;
  std::cout << std::endl;


  std::string specialCases( "" );
  std::string validBlocks( "" );
  std::string minimumScaleType( "" );
  std::string minimumScaleArgument( "" );
  std::string fixedScaleType( "" );
  std::string fixedScaleArgument( "" );
  std::string maximumScaleType( "" );
  std::string maximumScaleArgument( "" );
  while( lhaManagerParser.readNextElement() )
  {
    if( lhaManagerParser.currentElementNameMatches( "SpecialCases" ) )
    {
      specialCases = lhaManagerParser.getTrimmedCurrentElementContent();
    }
    else if( lhaManagerParser.currentElementNameMatches( "ValidBlocks" ) )
    {
      validBlocks = lhaManagerParser.getTrimmedCurrentElementContent();
    }
    else if( lhaManagerParser.currentElementNameMatches(
                                              "RenormalizationScaleChoices" ) )
    {
      BOL::AsciiXmlParser scaleChoiceParser;
      scaleChoiceParser.loadString(
                          lhaManagerParser.getTrimmedCurrentElementContent() );
      while( scaleChoiceParser.readNextElement() )
      {
        BOL::AsciiXmlParser elementParser;
        elementParser.loadString(
                         scaleChoiceParser.getTrimmedCurrentElementContent() );
        std::string evaluationType( "" );
        std::string evaluationArgument( "" );
        while( elementParser.readNextElement() )
        {
          if( elementParser.currentElementNameMatches( "EvaluationType" ) )
          {
            evaluationType = elementParser.getTrimmedCurrentElementContent();
          }
          else if( elementParser.currentElementNameMatches(
                                                       "EvaluationArgument" ) )
          {
            evaluationArgument
            = elementParser.getTrimmedCurrentElementContent();
          }
        }
        if( scaleChoiceParser.currentElementNameMatches(
                                                        "MinimumScaleBound" ) )
        {
          minimumScaleType = evaluationType;
          minimumScaleArgument = evaluationArgument;
        }
        else if( scaleChoiceParser.currentElementNameMatches(
                                                         "FixedScaleChoice" ) )
        {
          fixedScaleType = evaluationType;
          fixedScaleArgument = evaluationArgument;
        }
        else if( scaleChoiceParser.currentElementNameMatches(
                                                        "MaximumScaleBound" ) )
        {
          maximumScaleType = evaluationType;
          maximumScaleArgument = evaluationArgument;
        }
      }
    }
  }
  lhaManagerParser.closeFile();

  // Overwriting scale choices for RGE-improved comparison with fixed-case:
  minimumScaleType = fixedScaleType = maximumScaleType = "FixedNumber";
  minimumScaleArgument = fixedScaleArgument = maximumScaleArgument = "1000.0";

  VevaciousPlusPlus::LesHouchesAccordBlockEntryManager* lhaParameterManager(
                         new VevaciousPlusPlus::SlhaCompatibleWithSarahManager(
                                                                   validBlocks,
                                                              minimumScaleType,
                                                          minimumScaleArgument,
                                                                fixedScaleType,
                                                            fixedScaleArgument,
                                                              maximumScaleType,
                                                      maximumScaleArgument ) );

  LHPC::SlhaSimplisticInterpreter oldParser;
  oldParser.readFile( slhaFileName );
  LHPC::SimpleLhaParser newParser;
  newParser.ReadFile( slhaFileName );
  std::string inputString( "" );
  while( inputString != "q" )
  {
    std::cout << std::endl << "Input: ";
    std::getline( std::cin, inputString );
    std::list< std::pair< double, std::string > >
    oldValues( oldParser.getScalesPairedWithValues( inputString ) );
    std::list< std::pair< std::string, double > >
    newValues;
    newParser( inputString,
               newValues );
    std::cout << std::endl << "Old output: { ";
    for( std::list< std::pair< double, std::string > >::const_iterator
         oldPair( oldValues.begin() );
         oldPair != oldValues.end();
         ++oldPair )
    {
      if( oldPair != oldValues.begin() )
      {
        std::cout << "," << std::endl;
      }
      std::cout
      << "[ " << oldPair->first << ", \"" << oldPair->second << "\" ]";
    }
    std::cout << " }";
    std::cout << std::endl << "New output: { ";
    for( std::list< std::pair< std::string, double > >::const_iterator
         newPair( newValues.begin() );
         newPair != newValues.end();
         ++newPair )
    {
      if( newPair != newValues.begin() )
      {
        std::cout << "," << std::endl;
      }
      std::cout
      << "[ \"" << newPair->first << "\", " << newPair->second << " ]";
    }
  }



  // Cleaning up.
  delete lhaParameterManager;

  // debugging:
  /**/std::cout << std::endl << "debugging:"
  << std::endl
  << "End of testing old and new LHA parsers.";
  std::cout << std::endl;/**/

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

  std::cout
  << std::endl
  << "Vevacious finished running.";
  std::cout << std::endl;

  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}
