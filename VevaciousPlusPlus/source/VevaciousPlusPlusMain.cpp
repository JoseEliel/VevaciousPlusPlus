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

  // debugging:
  /**/std::cout << std::endl << "debugging:"
  << std::endl
  << "Testing old and new potential functions.";
  std::cout << std::endl;/**/

  std::string const slhaFileName( "CMSSM_CCB.slha.out" );

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
  while( lhaManagerParser.readNextElement() )
  {
    std::cout
    << std::endl
    << "lhaManagerParser.getCurrentElementName() = "
    << lhaManagerParser.getCurrentElementName();
    std::cout << std::endl;

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
        if( scaleChoiceParser.currentElementNameMatches( "FixedScaleChoice" ) )
        {
          fixedScaleType = evaluationType;
          fixedScaleArgument = evaluationArgument;
        }
        else if( scaleChoiceParser.currentElementNameMatches(
                                                        "MinimumScaleBound" ) )
        {
          minimumScaleType = evaluationType;
          minimumScaleArgument = evaluationArgument;
        }
      }
    }
  }
  lhaManagerParser.closeFile();

  std::cout
  << std::endl << "specialCases = \"" << specialCases << "\""
  << std::endl << "validBlocks = \"" << validBlocks << "\""
  << std::endl << "fixedScaleType = \"" << fixedScaleType << "\""
  << std::endl << "fixedScaleArgument = \"" << fixedScaleArgument << "\""
  << std::endl << "minimumScaleType = \"" << minimumScaleType << "\""
  << std::endl << "minimumScaleArgument = \"" << minimumScaleArgument << "\"";
  std::cout << std::endl;


  VevaciousPlusPlus::LesHouchesAccordBlockEntryManager*
  lhaParameterManager( NULL );
  if( specialCases == "SarahMssm" )
  {
    lhaParameterManager
    = new VevaciousPlusPlus::SlhaCompatibleWithSarahManager( validBlocks,
                                                             minimumScaleType,
                                                          minimumScaleArgument,
                                                             fixedScaleType,
                                                          fixedScaleArgument );
  }
  else if( specialCases == "SlhaMssm" )
  {
    lhaParameterManager
    = new VevaciousPlusPlus::SlhaBlocksWithSpecialCasesManager( validBlocks,
                                                              minimumScaleType,
                                                          minimumScaleArgument,
                                                                fixedScaleType,
                                                          fixedScaleArgument );
  }
  else
  {
    lhaParameterManager
    = new VevaciousPlusPlus::LesHouchesAccordBlockEntryManager( validBlocks,
                                                              minimumScaleType,
                                                          minimumScaleArgument,
                                                                fixedScaleType,
                                                          fixedScaleArgument );
  }

  std::string const oldModelFilename( "RealMssmWithStauAndStopVevs.vin" );

  VevaciousPlusPlus::RunningParameterManager
  slhaManager( *lhaParameterManager );
  VevaciousPlusPlus::OldFixedScaleOneLoopPotential
  oldFixedScale( oldModelFilename,
                 10.0,
                 true,
                 0.5,
                 slhaManager );

  VevaciousPlusPlus::OldRgeImprovedOneLoopPotential
  oldRgeImproved( oldFixedScale );

  slhaManager.UpdateSlhaData( slhaFileName );
  oldFixedScale.UpdateSelfForNewSlha( slhaManager );
  oldRgeImproved.UpdateSelfForNewSlha( slhaManager );

  std::string const
  newModelFilename( "NewFormatRealMssmWithStauAndStopVevsPotential.vin" );
  VevaciousPlusPlus::FixedScaleOneLoopPotential
  newFixedScale( newModelFilename,
                 0.5,
                 *lhaParameterManager );

  VevaciousPlusPlus::RgeImprovedOneLoopPotential
  newRgeImproved( newFixedScale );


  lhaParameterManager->NewParameterPoint( slhaFileName );

  std::cout
  << std::endl
  << "Parameter manager = " << lhaParameterManager->AsDebuggingString();
  std::cout << std::endl;


  double const oldFixedOriginTree( oldFixedScale.QuickApproximation(
                                         oldFixedScale.FieldValuesOrigin() ) );
  double const
  oldFixedOriginLoop( oldFixedScale( oldFixedScale.FieldValuesOrigin() ) );

  double const oldRgeOriginTree( oldRgeImproved.QuickApproximation(
                                        oldRgeImproved.FieldValuesOrigin() ) );
  double const
  oldRgeOriginLoop( oldRgeImproved( oldRgeImproved.FieldValuesOrigin() ) );

  double const newFixedOriginTree( newFixedScale.PolynomialApproximation()(
                                         newFixedScale.FieldValuesOrigin() ) );
  double const
  newFixedOriginLoop( newFixedScale( newFixedScale.FieldValuesOrigin() ) );

  double const newRgeOriginTree( newRgeImproved.PolynomialApproximation()(
                                        newRgeImproved.FieldValuesOrigin() ) );
  double const
  newRgeOriginLoop( newRgeImproved( newRgeImproved.FieldValuesOrigin() ) );

  std::cout  << std::endl
  << "oldFixedOriginTree = " << oldFixedOriginTree
  << ", newFixedOriginTree = " << newFixedOriginTree << std::endl
  << "oldRgeOriginTree = " << oldRgeOriginTree
  << ", newRgeOriginTree = " << newRgeOriginTree << std::endl
  << "oldFixedOriginLoop = " << oldFixedOriginLoop
  << ", newFixedOriginLoop = " << newFixedOriginLoop << std::endl
  << "oldRgeOriginLoop = " << oldRgeOriginLoop
  << ", newRgeOriginLoop = " << newRgeOriginLoop << std::endl;
  std::cout << std::endl;


  // Cleaning up.
  delete lhaParameterManager;

  // debugging:
  /**/std::cout << std::endl << "debugging:"
  << std::endl
  << "End of testing old and new potential functions.";
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
