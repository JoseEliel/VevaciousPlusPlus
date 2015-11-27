/*
 * VevaciousPlusPlusMain.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */


#include "CommonIncludes.hpp"
#include "VevaciousPlusPlus.hpp"
#include "LagrangianParameterManagement/RunningParameterManager.hpp"
#include "PotentialEvaluation/PotentialFunctions/OldFixedScaleOneLoopPotential.hpp"
#include "PotentialEvaluation/PotentialFunctions/OldRgeImprovedOneLoopPotential.hpp"
#include "PotentialMinimization/HomotopyContinuation/OldHom4ps2Runner.hpp"
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

  std::string const oldModelFilename( "RealMssmWithStauAndStopVevs.vin" );

  LHPC::SlhaSimplisticInterpreter oldParser;
  oldParser.readFile( slhaFileName );
  LHPC::SimpleLhaParser newParser;
  newParser.ReadFile( slhaFileName );
  std::string inputString( "" );
  while( inputString != "q" )
  {
    std::cout << std::endl << "Input: ";
    std::cin >> inputString;
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



  VevaciousPlusPlus::RunningParameterManager
  slhaManager( *lhaParameterManager );
  slhaManager.UpdateSlhaData( slhaFileName );
  VevaciousPlusPlus::OldFixedScaleOneLoopPotential
  oldFixedScale( oldModelFilename,
                 20.0,
                 true,
                 0.5,
                 slhaManager );
  std::vector< double > const fieldOrigin( oldFixedScale.FieldValuesOrigin() );

  VevaciousPlusPlus::OldRgeImprovedOneLoopPotential
  oldRgeImproved( oldFixedScale );

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
  std::vector< double > fixedParameterValues;
  lhaParameterManager->ParameterValues(
                     log( lhaParameterManager->AppropriateSingleFixedScale() ),
                                        fixedParameterValues );
  double const oldFixedOriginLoop( oldFixedScale( fieldOrigin ) );
  double const oldRgeOriginLoop( oldRgeImproved( fieldOrigin ) );
  double const newFixedOriginLoop( newFixedScale( fieldOrigin ) );
  double const newRgeOriginLoop( newRgeImproved( fieldOrigin ) );

  /*
  oldFixedScale.WriteAsPython( "OldFixedPython.py" );
  newFixedScale.WriteAsPython( "NewFixedPython.py" );
  */

  /*std::cout
  << std::endl
  << "Parameter manager = " << lhaParameterManager->AsDebuggingString();
  std::cout << std::endl;
  std::cout << "Old fixed tree = "
   << oldFixedScale.TreeLevelPotential().AsDebuggingString();
  std::cout << std::endl;
  std::cout << "New fixed tree = "
  << newFixedScale.PolynomialApproximation().AsDebuggingString();
  std::cout << std::endl;
  std::cout << "Old RGE tree = "
   << oldRgeImproved.TreeLevelPotential().AsDebuggingString();
  std::cout << std::endl;
  std::cout << "New RGE tree = "
  << newRgeImproved.PolynomialApproximation().AsDebuggingString();
  std::cout << std::endl;


  std::vector< double > unitHd( fieldOrigin );
  unitHd.front() = 1.0;
  std::vector< double > fixedParameterValues;
  lhaParameterManager->ParameterValues(
                     log( lhaParameterManager->AppropriateSingleFixedScale() ),
                                        fixedParameterValues );

  double const
  oldFixedOriginTree( oldFixedScale.QuickApproximation( fieldOrigin ) );
  double const
  oldRgeOriginTree( oldRgeImproved.QuickApproximation( fieldOrigin ) );
  double const
  newFixedOriginTree( newFixedScale.PolynomialApproximation()( fieldOrigin ) );
  double const newRgeOriginTree( newRgeImproved.PolynomialApproximation()(
                                                          fixedParameterValues,
                                                               fieldOrigin ) );

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

  std::cout  << std::endl
  << "oldFixed tree: at origin = "
  << oldFixedScale.QuickApproximation( fieldOrigin )
  << ", at unitHd = " << oldFixedScale.QuickApproximation( unitHd )
  << std::endl
  << "oldRge tree: at origin = "
  << oldRgeImproved.QuickApproximation( fieldOrigin )
  << ", at unitHd = " << oldRgeImproved.QuickApproximation( unitHd )
  << std::endl
  << "newFixed tree: at origin = "
  << newFixedScale.PolynomialApproximation()( fieldOrigin )
  << ", at unitHd = " << newFixedScale.PolynomialApproximation()( unitHd )
  << std::endl
  << "newRge tree: at origin = "
  << newRgeImproved.PolynomialApproximation()( fixedParameterValues,
                                               fieldOrigin )
  << ", at unitHd = "
  << newRgeImproved.PolynomialApproximation()( fixedParameterValues,
                                               unitHd )
  << std::endl;

  std::cout << std::endl
  << "Fixed scale parameter values = { ";
  for( std::vector< double >::const_iterator
       parameterValue( fixedParameterValues.begin() );
       parameterValue < fixedParameterValues.end();
       ++parameterValue )
  {
    std::cout << *parameterValue << ", ";
  }
  std::cout << "}"<< std::endl;
  std::cout
  << std::endl
  << "Random field configurations!" << std::endl;
  std::vector< double > randomConfiguration( fieldOrigin );
  for( unsigned int randomCount( 0 );
       randomCount < 10;
       ++randomCount )
  {
    for( std::vector< double >::iterator
         fieldValue( randomConfiguration.begin() );
         fieldValue < randomConfiguration.end();
         ++fieldValue )
    {
      *fieldValue
      = BOL::UsefulStuff::flatRandomDouble( ( -200.0 * randomCount ),
                                            ( 200.0 * randomCount ) );
    }

    std::cout
    << std::endl
    << oldFixedScale.FieldConfigurationAsMathematica( randomConfiguration )
    << std::endl
    << " - level\t\tOld fixed\t\tNew fixed\t\tOld RGE  \t\tNew RGE  "
    << std::endl
    << " - tree:\t\t"
    << oldFixedScale.QuickApproximation( randomConfiguration )
    << "\t\t"
    << newFixedScale.PolynomialApproximation()( randomConfiguration )
    << "\t\t"
    << oldRgeImproved.QuickApproximation( randomConfiguration )
    << "\t\t"
    << newRgeImproved.PolynomialApproximation()( fixedParameterValues,
                                                 randomConfiguration )
    << std::endl;
    std::cout
    << " - loop abs:\t\t"
    << oldFixedRandomLoop
    << "\t\t"
    << newFixedRandomLoop
    << "\t\t"
    << oldRgeRandomLoop
    << "\t\t"
    << newRgeRandomLoop;
    std::cout << std::endl;
    std::cout
    << " - loop rel:\t\t"
    << ( oldFixedRandomLoop - oldFixedOriginLoop )
    << "\t\t"
    << ( newFixedRandomLoop - newFixedOriginLoop )
    << "\t\t"
    << ( oldRgeRandomLoop - oldRgeOriginLoop )
    << "\t\t"
    << ( newRgeRandomLoop - newRgeOriginLoop );
    std::cout << std::endl;

  }
  std::cout << std::endl;
  std::cout
  << std::endl
  << "Random field configurations at random temperatures!" << std::endl;
  double randomTemperature( 0.0 );
  for( unsigned int randomCount( 0 );
       randomCount < 10;
       ++randomCount )
  {
    for( std::vector< double >::iterator
         fieldValue( randomConfiguration.begin() );
         fieldValue < randomConfiguration.end();
         ++fieldValue )
    {
      *fieldValue
      = BOL::UsefulStuff::flatRandomDouble( ( -200.0 * randomCount ),
                                            ( 200.0 * randomCount ) );
    }
    randomTemperature
    = BOL::UsefulStuff::flatRandomDouble( 0.0,
                                          ( 100.0 * randomCount ) );

    double const oldFixedRandomLoop( oldFixedScale( randomConfiguration,
                                                    randomTemperature ) );
    double const newFixedRandomLoop( newFixedScale( randomConfiguration,
                                                    randomTemperature ) );
    double const oldRgeRandomLoop( oldRgeImproved( randomConfiguration,
                                                   randomTemperature ) );
    double const newRgeRandomLoop( newRgeImproved( randomConfiguration,
                                                   randomTemperature ) );

    std::cout
    << std::endl
    << oldFixedScale.FieldConfigurationAsMathematica( randomConfiguration )
    << " T = " << randomTemperature << std::endl
    << " - level\t\tOld fixed\t\tNew fixed\t\tOld RGE  \t\tNew RGE  "
    << std::endl
    << " - tree:\t\t"
    << oldFixedScale.QuickApproximation( randomConfiguration )
    << "\t\t"
    << newFixedScale.PolynomialApproximation()( randomConfiguration )
    << "\t\t"
    << oldRgeImproved.QuickApproximation( randomConfiguration )
    << "\t\t"
    << newRgeImproved.PolynomialApproximation()( fixedParameterValues,
                                                 randomConfiguration )
    << std::endl;
    std::cout
    << " - loop abs:\t\t"
    << oldFixedRandomLoop
    << "\t\t"
    << newFixedRandomLoop
    << "\t\t"
    << oldRgeRandomLoop
    << "\t\t"
    << newRgeRandomLoop;
    std::cout << std::endl;
    std::cout
    << " - loop rel:\t\t"
    << ( oldFixedRandomLoop - oldFixedOriginLoop )
    << "\t\t"
    << ( newFixedRandomLoop - newFixedOriginLoop )
    << "\t\t"
    << ( oldRgeRandomLoop - oldRgeOriginLoop )
    << "\t\t"
    << ( newRgeRandomLoop - newRgeOriginLoop );
    std::cout << std::endl;

  }
  std::cout << std::endl;*/

  /*
  std::string const
  pathToHom4ps2( "/home/bol/BOL/ProjectDependencies/HOM4PS2/" );
  std::string const homotopyType( "2" );

  std::vector< std::vector< double > > oldHom4ps2Results;
  VevaciousPlusPlus::OldHom4ps2Runner oldHom4ps2Runner(
                         *(oldFixedScale.HomotopyContinuationTargetSystem()),
                                                        pathToHom4ps2,
                                                        homotopyType );
  oldFixedScale.HomotopyContinuationTargetSystem()->UpdateSelfForNewSlha(
                                                                 slhaManager );
  oldHom4ps2Runner( oldHom4ps2Results );

  std::vector< std::vector< double > > newHom4ps2Results;
  VevaciousPlusPlus::PolynomialSystemSolver*
  newHom4ps2Runner( new VevaciousPlusPlus::Hom4ps2Runner( pathToHom4ps2,
                                                          homotopyType,
                                                          1.0 ) );
  VevaciousPlusPlus::PolynomialAtFixedScalesSolver*
  newSolver( new VevaciousPlusPlus::PolynomialAtFixedScalesSolver(
                                      newFixedScale.PolynomialApproximation(),
                                                          *lhaParameterManager,
                                                              newHom4ps2Runner,
                                                                   1,
                                                                   true,
                                    newFixedScale.NumberOfFieldVariables() ) );
  (*newSolver)( newHom4ps2Results );

  std::list< std::vector< double > > vectorSorter( oldHom4ps2Results.begin(),
                                                   oldHom4ps2Results.end() );
  vectorSorter.sort();
  oldHom4ps2Results.assign( vectorSorter.begin(),
                            vectorSorter.end() );
  vectorSorter.assign( newHom4ps2Results.begin(),
                       newHom4ps2Results.end() );
  vectorSorter.sort();
  newHom4ps2Results.assign( vectorSorter.begin(),
                            vectorSorter.end() );

  size_t const numberOfSolutions( std::max( oldHom4ps2Results.size(),
                                            newHom4ps2Results.size() ) );
  size_t const numberOfFields( fieldOrigin.size() );

  for( size_t solutionIndex( 0 );
       solutionIndex < numberOfSolutions;
       ++solutionIndex )
  {
    std::cout
    << std::endl
    << "Old\t\t\tNew"
    << std::endl
    << "{ ";
    for( size_t fieldIndex( 0 );
         fieldIndex < newFixedScale.NumberOfFieldVariables();
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << std::endl;
        std::cout << "  ";
      }

      if( solutionIndex < oldHom4ps2Results.size() )
      {
        std::cout << oldHom4ps2Results[ solutionIndex ][ fieldIndex ];
      }
      else
      {
        std::cout << "missing";
      }
      if( fieldIndex < ( numberOfFields - 1 ) )
      {
        std::cout << ",";
      }
      std::cout << "\t\t";
      if( solutionIndex < newHom4ps2Results.size() )
      {
        std::cout << newHom4ps2Results[ solutionIndex ][ fieldIndex ];
      }
      else
      {
        std::cout << "missing";
      }
      if( fieldIndex < ( numberOfFields - 1 ) )
      {
        std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << " }" << std::endl;
  }

  VevaciousPlusPlus::MinuitPotentialMinimizer oldRoller( oldFixedScale,
                                                         0.1,
                                                         1.0,
                                                         1 );
  VevaciousPlusPlus::MinuitPotentialMinimizer*
  newRoller( new VevaciousPlusPlus::MinuitPotentialMinimizer( newFixedScale,
                                                              0.1,
                                                              1.0,
                                                              1 ) );

  for( std::vector< std::vector< double > >::const_iterator
       treeSolution( newHom4ps2Results.begin() );
       treeSolution != newHom4ps2Results.end();
       ++treeSolution )
  {
    VevaciousPlusPlus::PotentialMinimum
    oldMinimum( oldRoller( *treeSolution ) );
    VevaciousPlusPlus::PotentialMinimum
    newMinimum( (*newRoller)( *treeSolution ) );
    std::cout
    << std::endl
    << "Starting point \t\t\t\t"
    << newFixedScale.FieldConfigurationAsMathematica( *treeSolution )
    << std::endl << "New roller with old potential => \t"
    << newFixedScale.FieldConfigurationAsMathematica(
                                              oldMinimum.FieldConfiguration() )
    << std::endl << "New roller with new potential => \t"
    << newFixedScale.FieldConfigurationAsMathematica(
                                              newMinimum.FieldConfiguration() )
    << std::endl << "Potential difference from origin: old(old) = "
    << ( oldFixedScale( oldMinimum.FieldConfiguration() )
         - oldFixedOriginLoop )
    << ", old(new) = " << ( oldFixedScale( newMinimum.FieldConfiguration() )
                            - oldFixedOriginLoop )
    << ", new(old) = " << ( newFixedScale( oldMinimum.FieldConfiguration() )
                            - newFixedOriginLoop )
    << ", new(new) = " << ( newFixedScale( newMinimum.FieldConfiguration() )
                            - newFixedOriginLoop );
    std::cout << std::endl;
  }

  std::vector< double > fieldConfiguration( newHom4ps2Results.back() );
  std::cout
  << std::endl
  << oldFixedScale.PrintEvaluation( fieldConfiguration )
  << std::endl
  << newFixedScale.PrintEvaluation( fieldConfiguration );
  std::cout << std::endl;



  VevaciousPlusPlus::GradientFromStartingPoints newMinimizer( newFixedScale,
                                                              newSolver,
                                                              newRoller,
                                                              0.05,
                                                              4.0 );

  newMinimizer.FindMinima( 0.0 );
  */

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
