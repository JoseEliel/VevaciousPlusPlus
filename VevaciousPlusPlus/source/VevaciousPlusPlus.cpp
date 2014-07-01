/*
 * VevaciousPlusPlus.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  VevaciousPlusPlus::VevaciousPlusPlus(
                                  std::string const& initializationFileName ) :
    potentialFunction( NULL ),
    potentialFunctionDeleter( NULL ),
    runningParameterManager(),
    potentialMinimizer( NULL ),
    potentialMinimizerDeleter( NULL ),
    tunnelingCalculator( NULL ),
    tunnelingCalculatorDeleter( NULL ),
    currentTime()
  {
    BOL::AsciiXmlParser fileParser;
    BOL::AsciiXmlParser elementParser;
    fileParser.openRootElementOfFile( initializationFileName );
    while( fileParser.readNextElement() )
    {
      elementParser.loadString( fileParser.getTrimmedCurrentElementContent() );
      if( fileParser.currentElementNameMatches( "PotentialClass" ) )
      {
        // <PotentialClass> should have child elements <ClassType> and
        // <ConstructorArguments>.
        elementParser.readNextElement();
        std::string
        potentialClass( elementParser.getTrimmedCurrentElementContent() );
        elementParser.readNextElement();
        std::string constructorArguments(
                             elementParser.getTrimmedCurrentElementContent() );
        if( potentialClass.compare( "FixedScaleOneLoopPotential" ) == 0 )
        {
          potentialFunction
          = new FixedScaleOneLoopPotential( constructorArguments,
                                            runningParameterManager );
        }
        else if( potentialClass.compare( "RgeImprovedOneLoopPotential" ) == 0 )
        {
          potentialFunction
          = new RgeImprovedOneLoopPotential( constructorArguments,
                                             runningParameterManager );
        }
        else
        {
          std::stringstream errorStream;
          errorStream
          << "<PotentialClass> was not a recognized form! The only types"
          << " currently valid are \"FixedScaleOneLoopPotential\" and"
          << " \"RgeImprovedOneLoopPotential\".";
          throw std::runtime_error( errorStream.str() );
        }
        potentialFunctionDeleter = potentialFunction;
      }
      else if( fileParser.currentElementNameMatches( "MinimizerClass" ) )
      {
        // <MinimizerClass> should have child elements <ClassName> and
        // <ConstructorArguments>.
        elementParser.readNextElement();
        std::string minimizerClass(
                             elementParser.getTrimmedCurrentElementContent() );
        elementParser.readNextElement();
        std::string constructorArguments(
                             elementParser.getTrimmedCurrentElementContent() );
        if( minimizerClass.compare( "GradientFromStartingPoints" ) == 0 )
        {
          potentialMinimizer
          = new GradientFromStartingPoints( potentialFunction,
                                            constructorArguments,
                                            runningParameterManager );
        }
        else
        {
          std::stringstream errorStream;
          errorStream
          << "<MinimizerClass> was not a recognized form! The only type"
          << " currently valid is \"GradientFromStartingPoints\".";
          throw std::runtime_error( errorStream.str() );
        }
        potentialMinimizerDeleter = potentialMinimizer;
      }
    }


    VevaciousPlusPlus::HomotopyContinuationSolver*
    homotopyContinuationSolver( NULL );

    std::string homotopyContinuationClass(
                             argumentParser.fromTag( "HomotopyContinuationClass",
                                       "BasicPolynomialHomotopyContinuation" ) );
    if( homotopyContinuationClass.compare(
                                   "BasicPolynomialHomotopyContinuation" ) == 0 )
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
                                                               "2" ) );

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
      std::string numberOfNodesString(
                       argumentParser.fromTag( "NumberOfMinuitNodesForTunneling",
                                               "10" ) );
      double numberOfNodesDouble( -1.0 );
      validInput = BOL::StringParser::stringIsDouble( numberOfNodesString,
                                                      numberOfNodesDouble );
      if( !validInput
          ||
          !( numberOfNodesDouble >= 1.0 ) )
      {
        std::cout
        << std::endl
        << "NumberOfMinuitNodesForTunneling was not a number > 1, but was "
        << numberOfNodesString << "."
        << " Aborting!";
        std::cout << std::endl;

        return EXIT_FAILURE;
      }

      std::string initialStepSizeString(
             argumentParser.fromTag( "InitialMinuitStepSizeFractionForTunneling",
                                     "0.1" ) );
      double initialStepSize( -1.0 );
      validInput = BOL::StringParser::stringIsDouble( initialStepSizeString,
                                                      initialStepSize );
      if( !validInput
          ||
          !( initialStepSize > 0.0 ) )
      {
        std::cout
        << std::endl
        << "NumberOfMinuitNodesForTunneling was not a number > 0, but was "
        << initialStepSizeString << "."
        << " Aborting!";
        std::cout << std::endl;

        return EXIT_FAILURE;
      }

      tunnelingCalculator
      = new VevaciousPlusPlus::MinuitBounceActionMinimizer( *potentialFunction,
                                                            tunnelingStrategy,
                                                    survivalProbabilityThreshold,
                                                     (size_t)numberOfNodesDouble,
                                                            initialStepSize );

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

  }

  VevaciousPlusPlus::VevaciousPlusPlus( BOL::ArgumentParser& argumentParser,
                                        SlhaManager& slhaManager,
                                        PotentialMinimizer& potentialMinimizer,
                                   TunnelingCalculator& tunnelingCalculator ) :
    potentialFunction( NULL ),
    potentialFunctionDeleter( NULL ),
    runningParameterManager(),
    potentialMinimizer( &potentialMinimizer ),
    potentialMinimizerDeleter( NULL ),
    tunnelingCalculator( &tunnelingCalculator ),
    tunnelingCalculatorDeleter( NULL ),
    currentTime()
  {
    // This constructor is just an initialization list.
  }

  VevaciousPlusPlus::~VevaciousPlusPlus()
  {
    delete tunnelingCalculatorDeleter;
    delete potentialMinimizerDeleter;
    delete potentialFunctionDeleter;
  }


  void VevaciousPlusPlus::RunPoint( std::string const& parameterFilename )
  {
    time( &currentTime );
    std::cout
    << std::endl
    << "Running " << parameterFilename << " starting at "
    << ctime( &currentTime );
    std::cout << std::endl;
    BOL::BasicTimer stageTimer;
    BOL::BasicTimer totalTimer;

    runningParameterManager.UpdateSlhaData( parameterFilename );
    potentialMinimizer->FindMinima( 0.0 );

    time( &currentTime );
    std::cout
    << std::endl
    << "Minimization of potential took " << stageTimer.secondsSinceStart()
    << " seconds, finished at " << ctime( &currentTime );
    std::cout << std::endl;

    if( potentialMinimizer->DsbVacuumIsMetastable() )
    {
      stageTimer.setStartTime();
      tunnelingCalculator->CalculateTunneling( potentialMinimizer->DsbVacuum(),
                                           potentialMinimizer->PanicVacuum() );

      time( &currentTime );
      std::cout
      << std::endl
      << "Tunneling calculation took " << stageTimer.secondsSinceStart()
      << " seconds, finished at " << ctime( &currentTime );
      std::cout << std::endl;
      std::cout << std::endl;
    }

    time( &currentTime );
    std::cout
    << std::endl
    << "Total running time was " << totalTimer.secondsSinceStart()
    << " seconds, finished at " << ctime( &currentTime );
    std::cout << std::endl;
  }

  void VevaciousPlusPlus::WriteXmlResults( std::string const& xmlFilename )
  {
    std::ofstream xmlFile( xmlFilename.c_str() );
    xmlFile << "<VevaciousResults>\n"
    "  <ReferenceData>\n"
    "     <VevaciousVersion>\n"
    "       " << versionString << "\n"
    "     </VevaciousVersion>\n"
    "     <CitationArticle>\n"
    "       " << citationString << "\n"
    "     </CitationArticle>\n"
    "     <ResultTimestamp>\n"
    "       "
    << std::string( ctime( &currentTime ) )
    << "     </ResultTimestamp>\n"
    "  </ReferenceData>\n"
    "  <StableOrMetastable>\n"
    "    ";
    if( potentialMinimizer->DsbVacuumIsMetastable() )
    {
      xmlFile << "meta";
    }
    xmlFile << "stable\n"
    "  </StableOrMetastable>\n"
    << potentialMinimizer->DsbVacuum().AsVevaciousXmlElement( "DsbVacuum",
                                            potentialMinimizer->FieldNames() );
    if( potentialMinimizer->DsbVacuumIsMetastable() )
    {
      xmlFile
      << potentialMinimizer->PanicVacuum().AsVevaciousXmlElement(
                                                                "PanicVacuum",
                                            potentialMinimizer->FieldNames() );
      if( tunnelingCalculator->QuantumSurvivalProbability() >= 0.0 )
      {
        xmlFile << "  <ZeroTemperatureDsbSurvival>\n"
        "    <DsbSurvivalProbability>\n"
        "      " << tunnelingCalculator->QuantumSurvivalProbability() << "\n"
        "    </DsbSurvivalProbability>\n"
        "    <LogOfMinusLogOfDsbSurvival>\n"
        "      " << tunnelingCalculator->LogOfMinusLogOfQuantumProbability()
        << " <!-- this = ln(-ln(P)), so P = e^(-e^this)) -->\n"
        "    </LogOfMinusLogOfDsbSurvival>\n"
        "    <DsbLifetime>\n"
        "      " << tunnelingCalculator->QuantumLifetimeInSeconds()
        << " <!-- in seconds; age of observed Universe is 4.3E+17s -->\n"
        << "    </DsbLifetime>\n"
        "  </ZeroTemperatureDsbSurvival>\n";
      }
      else
      {
        xmlFile << "  <!-- Survival probability at zero temperature not"
                                                          " calculated. -->\n";
      }
      if( tunnelingCalculator->ThermalSurvivalProbability() >= 0.0 )
      {
        xmlFile << "  <NonZeroTemperatureDsbSurvival>\n"
        "    <DsbSurvivalProbability>\n"
        "      " << tunnelingCalculator->ThermalSurvivalProbability() << "\n"
        "    </DsbSurvivalProbability>\n"
        "    <LogOfMinusLogOfDsbSurvival>\n"
        "      " << tunnelingCalculator->LogOfMinusLogOfThermalProbability()
        << " <!-- this = ln(-ln(P)), so P = e^(-e^this)) --> \n"
        "    </LogOfMinusLogOfDsbSurvival>\n"
        "    <DominantTunnelingTemperature>\n"
        "      "
        << tunnelingCalculator->DominantTemperatureInGigaElectronVolts()
        << " <!-- in GeV -->\n"
        << "    </DominantTunnelingTemperature>\n"
        "  </NonZeroTemperatureDsbSurvival>\n";
      }
      else
      {
        xmlFile << "  <!-- Survival probability at non-zero temperatures not"
                                                          " calculated. -->\n";
      }
    }
    xmlFile << "</VevaciousResults>\n";
    xmlFile.close();
  }

  void VevaciousPlusPlus::WriteSlhaResults( std::string const& slhaFilename,
                                            bool const writeWarnings )
  {
    BOL::StringParser const slhaIndexMaker( 3,
                                            ' ',
                                            9,
                                            3,
                                            "" );
    BOL::StringParser const slhaDoubleMaker( 9,
                                             ' ',
                                             9,
                                             3 );
    std::fstream outputFile( slhaFilename.c_str() );
    if( !(outputFile.good()) )
    {
      throw std::runtime_error( "Could not open \"" + slhaFilename
                                + "\" to append results." );
    }
    long endPosition( outputFile.seekg( 0,
                                        std::ios::end ).tellg() );
    while( '\n' == (char)(outputFile.seekg( (--endPosition) ).peek()) )
    {
      // this loop just brings the get pointer back to the char before the
      // last '\n'.
    }
    outputFile.seekp( (++endPosition) );
    // the put pointer is now about to overwrite the 1st '\n' of the sequence
    // of '\n' characters ending the file.
    outputFile << "\n"
    "BLOCK VEVACIOUSSTABILITY # Results from VevaciousPlusPlus\n"
    "# version " << versionString << ", documented in " << citationString
    << "\n"
    "# Results written " << std::string( ctime( &currentTime ) )
    << "# [index] [verdict int]\n"
    "  1  ";
    if( potentialMinimizer->DsbVacuumIsStable() )
    {
      outputFile << "1  # Stable DSB vacuum\n";
    }
    else
    {
      outputFile << "0  # Metastable DSB vacuum\n";
    }
    outputFile << "BLOCK VEVACIOUSZEROTEMPERATURE # Results at T = 0\n"
    "# [index] [verdict float]\n";
    if( tunnelingCalculator->QuantumSurvivalProbability() >= 0.0 )
    {
      outputFile <<  "  1  " << slhaDoubleMaker.doubleToString(
                            tunnelingCalculator->QuantumSurvivalProbability() )
      << "  # Probability of DSB vacuum surviving 4.3E17 seconds\n";
      outputFile << "  2  " << slhaDoubleMaker.doubleToString(
                              tunnelingCalculator->QuantumLifetimeInSeconds() )
      << "  # Tunneling time out of DSB vacuum in seconds\n"
      "  3  " << slhaDoubleMaker.doubleToString(
                     tunnelingCalculator->LogOfMinusLogOfQuantumProbability() )
      << "  # L = ln(-ln(P)), => P = e^(-e^L)\n";
    }
    else
    {
      outputFile << "  1  " << slhaDoubleMaker.doubleToString(
          tunnelingCalculator->QuantumSurvivalProbability() )
      << "  # Not calculated: ignore this number\n"
      "  2  " << slhaDoubleMaker.doubleToString(
                              tunnelingCalculator->QuantumLifetimeInSeconds() )
      << "  # Not calculated: ignore this number\n"
      "  3  " << slhaDoubleMaker.doubleToString(
                     tunnelingCalculator->LogOfMinusLogOfQuantumProbability() )
      << "  # Not calculated: ignore this number\n";
    }
    outputFile << "BLOCK VEVACIOUSNONZEROTEMPERATURE # Results at T != 0\n"
    "# [index] [verdict float]\n";
    if( tunnelingCalculator->ThermalSurvivalProbability() >= 0.0 )
    {
      outputFile <<  "  1  " << slhaDoubleMaker.doubleToString(
                            tunnelingCalculator->ThermalSurvivalProbability() )
      << "  # Probability of DSB vacuum surviving thermal tunneling\n";
      outputFile << "  2  " << slhaDoubleMaker.doubleToString(
                tunnelingCalculator->DominantTemperatureInGigaElectronVolts() )
      << "  # Dominant tunneling temperature in GeV\n"
      "  3  " << slhaDoubleMaker.doubleToString(
                     tunnelingCalculator->LogOfMinusLogOfThermalProbability() )
      << "  # L = ln(-ln(P)), => P = e^(-e^L)\n";
    }
    else
    {
      outputFile << "  1  " << slhaDoubleMaker.doubleToString(
                            tunnelingCalculator->ThermalSurvivalProbability() )
      << "  # Not calculated: ignore this number\n"
      "  2  " << slhaDoubleMaker.doubleToString(
                 tunnelingCalculator->DominantTemperatureInGigaElectronVolts() )
      << "  # Not calculated: ignore this number\n"
      "  3  " << slhaDoubleMaker.doubleToString(
                     tunnelingCalculator->LogOfMinusLogOfThermalProbability() )
      << "  # Not calculated: ignore this number\n";
    }
    outputFile
    << "BLOCK VEVACIOUSFIELDNAMES # Field names for each index\n"
    "# [index] [field name in \"\"]\n";
    std::vector< std::string > const&
    fieldNames( potentialMinimizer->FieldNames() );
    for( size_t fieldIndex( 0 );
         fieldIndex < fieldNames.size();
         ++fieldIndex )
    {
      outputFile << slhaIndexMaker.intToString( fieldIndex ) << "  \""
      << fieldNames[ fieldIndex ] << "\"\n";
    }
    outputFile << "BLOCK VEVACIOUSDSBVACUUM # VEVs for DSB vacuum in GeV\n"
    "# [index] [field VEV in GeV]\n";
    std::vector< double > const&
    dsbFields( potentialMinimizer->DsbVacuum().FieldConfiguration() );
    for( size_t fieldIndex( 0 );
         fieldIndex < fieldNames.size();
         ++fieldIndex )
    {
      outputFile << slhaIndexMaker.intToString( fieldIndex ) << "  "
      << slhaDoubleMaker.doubleToString( dsbFields[ fieldIndex ] ) << "  # "
      << fieldNames[ fieldIndex ] << "\n";
    }
    outputFile
    << "BLOCK VEVACIOUSPANICVACUUM # ";
    if( potentialMinimizer->DsbVacuumIsMetastable() )
    {
      outputFile << "VEVs for panic vacuum in GeV\n";
    }
    else
    {
      outputFile << "Stable DSB vacuum => repeating DSB VEVs\n";
    }
    outputFile << "# [index] [field VEV in GeV]\n";
    std::vector< double > const*
    panicFields( &(potentialMinimizer->DsbVacuum().FieldConfiguration()) );
    if( potentialMinimizer->DsbVacuumIsMetastable() )
    {
      panicFields = &(potentialMinimizer->PanicVacuum().FieldConfiguration());
    }
    for( size_t fieldIndex( 0 );
         fieldIndex < fieldNames.size();
         ++fieldIndex )
    {
      outputFile << slhaIndexMaker.intToString( fieldIndex ) << "  "
      << slhaDoubleMaker.doubleToString( (*panicFields)[ fieldIndex ] )
      << "  # " << fieldNames[ fieldIndex ] << "\n";
    }
    outputFile.close();
  }

} /* namespace VevaciousPlusPlus */
