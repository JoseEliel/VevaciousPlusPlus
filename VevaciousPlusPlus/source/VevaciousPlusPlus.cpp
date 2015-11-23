/*
 * VevaciousPlusPlus.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{
  // This is the constructor for those who know what they are doing, to allow
  // the main function of the program to decide the components and pass them
  // in to the constructor, allowing for custom components without having to
  // edit the VevaciousPlusPlus files. Those wishing to just use the default
  // possibilities can use the other constructor, which takes the name of an
  // initialization file, and then creates the components based on the data
  // in that file. Since this constructor leaves ownedPotentialFunction,
  // ownedPotentialMinimizer, and ownedPotentialMinimizer all as
  // NULL, and no other function sets them, there should be no problem with
  // the destructor calling delete on these pointers, as they do not get set to
  // point at the addresses of the given components.
  VevaciousPlusPlus::VevaciousPlusPlus( PotentialMinimizer& potentialMinimizer,
                                   TunnelingCalculator& tunnelingCalculator ) :
    lagrangianParameterManager( &(potentialMinimizer.GetPotentialFunction(
                                          ).GetLagrangianParameterManager()) ),
    ownedLagrangianParameterManager( NULL ),
    slhaManager( NULL ),
    // potentialFunction( &(potentialMinimizer.GetPotentialFunction()) ),
    ownedPotentialFunction( NULL ),
    oldPotentialFunction( NULL ),
    potentialMinimizer( &potentialMinimizer ),
    ownedPotentialMinimizer( NULL ),
    tunnelingCalculator( &tunnelingCalculator ),
    ownedTunnelingCalculator( NULL ),
    currentTime()
  {
    // This constructor is just an initialization list.
  }

  // This is the constructor that we expect to be used in normal use: it reads
  // in an initialization file in XML with name given by initializationFileName
  // and then assembles the appropriate components for the objects. The body of
  // the constructor is just a lot of statements reading in XML elements and
  // creating new instances of components. It'd be great to be able to use
  // std::unique_ptrs but we're sticking to allowing non-C++11-compliant
  // compilers.
  VevaciousPlusPlus::VevaciousPlusPlus(
                                  std::string const& initializationFileName ) :
    lagrangianParameterManager( NULL ),
    ownedLagrangianParameterManager( NULL ),
    slhaManager( NULL ),
    // potentialFunction( NULL ),
    ownedPotentialFunction( NULL ),
    potentialMinimizer( NULL ),
    ownedPotentialMinimizer( NULL ),
    tunnelingCalculator( NULL ),
    ownedTunnelingCalculator( NULL ),
    currentTime()
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Really have to sort out this constructor. Remember to remove"
    << " slhaManager!";
    std::cout << std::endl;/**/

    std::string potentialFunctionInitializationFilename( "error" );
    std::string potentialMinimizerInitializationFilename( "error" );
    std::string tunnelingCalculatorInitializationFilename( "error" );
    BOL::AsciiXmlParser fileParser;
    fileParser.openRootElementOfFile( initializationFileName );
    while( fileParser.readNextElement() )
    {
      if( fileParser.currentElementNameMatches(
                                      "PotentialFunctionInitializationFile" ) )
      {
        potentialFunctionInitializationFilename
        = fileParser.getTrimmedCurrentElementContent();
      }
      else if( fileParser.currentElementNameMatches(
                                     "PotentialMinimizerInitializationFile" ) )
      {
        potentialMinimizerInitializationFilename
        = fileParser.getTrimmedCurrentElementContent();
      }
      else if( fileParser.currentElementNameMatches(
                                    "TunnelingCalculatorInitializationFile" ) )
      {
        tunnelingCalculatorInitializationFilename
        = fileParser.getTrimmedCurrentElementContent();
      }
    }

    // Now we know the overall picture of what components to set up.
    // ALL SUB-COMPONENTS FOR potentialMinimizer AND tunnelingCalculator ARE
    // MEMORY-MANAGED BY THE COMPONENTS THEMSELVES! The VevaciousPlusPlus
    // destructor only deletes ownedTunnelingCalculator,
    // ownedPotentialMinimizer, ownedPotentialFunction, and
    // ownedLagrangianParameterManager, while this constructor is allocating
    // much more memory than that. Everything would be clear if we restricted
    // ourselves to requiring a C++11-compliant compiler, as then we could have
    // the constructors use std::unique_ptrs, but, alas, we're sticking to
    // C++98.
    std::pair< LesHouchesAccordBlockEntryManager*,
               PotentialFromPolynomialWithMasses* >
    lagrangianParameterManagerAndPotentialFunction(
                          CreateLagrangianParameterManagerAndPotentialFunction(
                                   potentialFunctionInitializationFilename ) );
    lagrangianParameterManager = ownedLagrangianParameterManager
    = lagrangianParameterManagerAndPotentialFunction.first;
    // potentialFunction =
    ownedPotentialFunction
    = lagrangianParameterManagerAndPotentialFunction.second;
    ownedPotentialMinimizer
    = CreatePotentialMinimizer( *ownedPotentialFunction,
                                potentialMinimizerInitializationFilename );
    ownedTunnelingCalculator
    = CreateTunnelingCalculator( tunnelingCalculatorInitializationFilename );
  }

  VevaciousPlusPlus::~VevaciousPlusPlus()
  {
    delete ownedTunnelingCalculator;
    delete ownedPotentialMinimizer;
    delete ownedPotentialFunction;
    delete oldPotentialFunction;
    delete slhaManager;
    delete ownedLagrangianParameterManager;
  }


  // This runs the point parameterized by newInput, which for the default
  // case gives the name of a file with the input parameters, but could in
  // principle itself contain all the necessary parameters.
  void VevaciousPlusPlus::RunPoint( std::string const& newInput )
  {
    time( &currentTime );
    std::cout
    << std::endl
    << "Running \"" << newInput << "\" starting at "
    << ctime( &currentTime );
    std::cout << std::endl;
    BOL::BasicTimer stageTimer;
    BOL::BasicTimer totalTimer;

    lagrangianParameterManager->NewParameterPoint( newInput );
    slhaManager->UpdateSlhaData( newInput );
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
      tunnelingCalculator->CalculateTunneling(
                                    potentialMinimizer->GetPotentialFunction(),
                                               potentialMinimizer->DsbVacuum(),
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

  // This writes the results as an XML file.
  void VevaciousPlusPlus::WriteXmlResults( std::string const& xmlFilename )
  {
    std::ofstream xmlFile( xmlFilename.c_str() );
    xmlFile << "<VevaciousResults>\n"
    "  <ReferenceData>\n"
    "     <VevaciousVersion>\n"
    "       " << VersionInformation::currentVersion << "\n"
    "     </VevaciousVersion>\n"
    "     <CitationArticle>\n"
    "       " << VersionInformation::currentCitation << "\n"
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
    std::vector< std::string > const&
    fieldNames( potentialMinimizer->GetPotentialFunction().FieldNames() );
    xmlFile << "stable\n"
    "  </StableOrMetastable>\n"
    << potentialMinimizer->DsbVacuum().AsVevaciousXmlElement( "DsbVacuum",
                                                              fieldNames );
    if( potentialMinimizer->DsbVacuumIsMetastable() )
    {
      xmlFile
      << potentialMinimizer->PanicVacuum().AsVevaciousXmlElement(
                                                                "PanicVacuum",
                                                                fieldNames );
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

    std::cout << std::endl << "Wrote results in XML in file \"" << xmlFilename
    << "\"." << std::endl;
  }

  // This writes the results as an SLHA file.
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
                                             3,
                                             "" );
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
    "# version " << VersionInformation::currentVersion << ", documented in "
    << VersionInformation::currentCitation
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
    fieldNames( potentialMinimizer->GetPotentialFunction().FieldNames() );
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
    std::cout << std::endl << "Wrote results in SLHA format at end of file \""
    << slhaFilename << "\"." << std::endl;
  }

  // This puts the content of the current element of xmlParser into
  // contentDestination, interpreted as a bool represented by
  // case-insensitive "yes/no" or "y/n" or "true/false" or "t/f" or "0/1",
  // if the element's name matches elementName. If the element content
  // doesn't match any valid input, contentDestination is left untouched.
  void VevaciousPlusPlus::InterpretElementIfNameMatches(
                                          BOL::AsciiXmlParser const& xmlParser,
                                                std::string const& elementName,
                                                     bool& contentDestination )
  {
    if( xmlParser.currentElementNameMatches( elementName ) )
    {
      std::string contentString( xmlParser.getTrimmedCurrentElementContent() );
      BOL::StringParser::transformToLowercase( contentString );
      if( ( contentString.compare( "yes" ) == 0 )
          ||
          ( contentString.compare( "y" ) == 0 )
          ||
          ( contentString.compare( "true" ) == 0 )
          ||
          ( contentString.compare( "t" ) == 0 )
          ||
          ( contentString.compare( "1" ) == 0 ) )
      {
        contentDestination = true;
      }
      else if( ( contentString.compare( "no" ) == 0 )
               ||
               ( contentString.compare( "n" ) == 0 )
               ||
               ( contentString.compare( "false" ) == 0 )
               ||
               ( contentString.compare( "f" ) == 0 )
               ||
               ( contentString.compare( "0" ) == 0 ) )
      {
        contentDestination = false;
      }
    }
  }

  // This interprets the given string as the appropriate element of the
  // TunnelingCalculator::TunnelingStrategy enum.
  TunnelingCalculator::TunnelingStrategy
  VevaciousPlusPlus::InterpretTunnelingStrategy(
                                               std::string& tunnelingStrategy )
  {
    BOL::StringParser::transformToLowercase( tunnelingStrategy );
    if( ( tunnelingStrategy.compare( "defaulttunneling" ) == 0 )
        ||
        ( tunnelingStrategy.compare( "thermalthenquantum" ) == 0 ) )
    {
      return TunnelingCalculator::ThermalThenQuantum;
    }
    else if( tunnelingStrategy.compare( "quantumthenthermal" ) == 0 )
    {
      return TunnelingCalculator::QuantumThenThermal;
    }
    else if( tunnelingStrategy.compare( "justthermal" ) == 0 )
    {
      return TunnelingCalculator::JustThermal;
    }
    else if( tunnelingStrategy.compare( "justquantum" ) == 0 )
    {
      return TunnelingCalculator::JustQuantum;
    }
    return TunnelingCalculator::NoTunneling;
  }

  // This decides on the derived class (based on GradientMinimizer) to use
  // to minimize the potential and constructs it with the arguments parsed
  // from constructorArguments.
  PotentialMinimizer* VevaciousPlusPlus::SetUpGradientFromStartingPoints(
                                      std::string const& constructorArguments )
  {
    // We need to assemble the components for a GradientFromStartingPoints
    // object: a StartingPointFinder and a GradientMinimizer. We need to find
    // out what derived classes to actually use.
    BOL::AsciiXmlParser elementParser;
    elementParser.loadString( constructorArguments );
    StartingPointFinder* startingPointFinder( NULL );
    std::string startingPointFinderClass( "Hom4ps2Runner" );
    std::string startingPointFinderArguments( "" );
    GradientMinimizer* gradientMinimizer( NULL );
    std::string gradientMinimizerClass( "MinuitPotentialMinimizer" );
    std::string gradientMinimizerArguments( "" );
    double extremumSeparationThresholdFraction( 0.05 );
    double nonDsbRollingToDsbScalingFactor( 4.0 );
    // The <ConstructorArguments> for this class should have child elements
    // <StartingPointFinderClass> and <GradientMinimizerClass>, and
    // optionally <ExtremumSeparationThresholdFraction> and
    // <NonDsbRollingToDsbScalingFactor>.
    while( elementParser.readNextElement() )
    {
      ReadClassAndArguments( elementParser,
                             "StartingPointFinderClass",
                             startingPointFinderClass,
                             startingPointFinderArguments );
      ReadClassAndArguments( elementParser,
                             "GradientMinimizerClass",
                             gradientMinimizerClass,
                             gradientMinimizerArguments );
      InterpretElementIfNameMatches( elementParser,
                                     "ExtremumSeparationThresholdFraction",
                                     extremumSeparationThresholdFraction );
      InterpretElementIfNameMatches( elementParser,
                                     "NonDsbRollingToDsbScalingFactor",
                                     nonDsbRollingToDsbScalingFactor );
    }
    if( startingPointFinderClass.compare( "Hom4ps2Runner" ) == 0 )
    {
      startingPointFinder = SetUpHom4ps2Runner( startingPointFinderArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<StartingPointFinderClass> was not a recognized form! The only"
      << " type currently valid is \"Hom4ps2Runner\".";
      throw std::runtime_error( errorStream.str() );
    }
    if( gradientMinimizerClass.compare( "MinuitPotentialMinimizer" ) == 0 )
    {
      gradientMinimizer
      = SetUpMinuitPotentialMinimizer( gradientMinimizerArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<GradientMinimizerClass> was not a recognized form! The only type"
      << " currently valid is \"MinuitPotentialMinimizer\".";
      throw std::runtime_error( errorStream.str() );
    }
    // Now we have the components for potentialMinimizer:
    ownedPotentialMinimizer =
    new GradientFromStartingPoints( *ownedPotentialFunction,
                                    startingPointFinder,
                                    gradientMinimizer,
                                    extremumSeparationThresholdFraction,
                                    nonDsbRollingToDsbScalingFactor );
    return ownedPotentialMinimizer;
  }

  // This parses arguments from constructorArguments and uses them to
  // construct a MinuitPotentialMinimizer instance to use to minimize the
  // potential from a set of starting points.
  GradientMinimizer* VevaciousPlusPlus::SetUpMinuitPotentialMinimizer(
                                      std::string const& constructorArguments )
  {
    double errorFraction( 0.1 );
    double errorMinimum( 1.0 );
    int minuitStrategy( 1 );
    BOL::AsciiXmlParser xmlParser;
    xmlParser.loadString( constructorArguments );
    while( xmlParser.readNextElement() )
    {
      InterpretElementIfNameMatches( xmlParser,
                                     "InitialStepSizeFraction",
                                     errorFraction );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinimumInitialStepSize",
                                     errorMinimum );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinuitStrategy",
                                     minuitStrategy );
    }
    return new MinuitPotentialMinimizer( *ownedPotentialFunction,
                                         errorFraction,
                                         errorMinimum,
                                         minuitStrategy );
  }

  // This parses arguments from constructorArguments and uses them to
  // construct a CosmoTransitionsRunner instance to use to calculate the
  // tunneling from the DSB vacuum to the panic vacuum, if possible.
  TunnelingCalculator* VevaciousPlusPlus::SetUpCosmoTransitionsRunner(
                                      std::string const& constructorArguments )
  {
    // The <ConstructorArguments> for this class should have child elements
    // <TunnelingStrategy>, <SurvivalProbabilityThreshold>,
    // <CriticalTemperatureAccuracy>, <EvaporationBarrierResolution>,
    // <PathToCosmotransitions>, <PathResolution>, <MaxInnerLoops>, and
    // <MaxOuterLoops>.
    std::string tunnelingStrategy( "ThermalThenQuantum" );
    double survivalProbabilityThreshold( 0.1 );
    int thermalStraightPathFitResolution( 5 );
    int temperatureAccuracy( 7 );
    std::string pathToCosmotransitions( "./cosmoTransitions/" );
    int resolutionOfDsbVacuum( 20 );
    double vacuumSeparationFraction( 0.2 );
    int maxInnerLoops( 10 );
    int maxOuterLoops( 10 );
    BOL::AsciiXmlParser xmlParser;
    xmlParser.loadString( constructorArguments );
    while( xmlParser.readNextElement() )
    {
      InterpretElementIfNameMatches( xmlParser,
                                     "TunnelingStrategy",
                                     tunnelingStrategy );
      InterpretElementIfNameMatches( xmlParser,
                                     "SurvivalProbabilityThreshold",
                                     survivalProbabilityThreshold );
      InterpretElementIfNameMatches( xmlParser,
                                     "ThermalActionResolution",
                                     thermalStraightPathFitResolution );
      InterpretElementIfNameMatches( xmlParser,
                                     "CriticalTemperatureAccuracy",
                                     temperatureAccuracy );
      InterpretElementIfNameMatches( xmlParser,
                                     "PathToCosmotransitions",
                                     pathToCosmotransitions );
      InterpretElementIfNameMatches( xmlParser,
                                     "PathResolution",
                                     resolutionOfDsbVacuum );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinimumVacuumSeparationFraction",
                                     vacuumSeparationFraction );
      InterpretElementIfNameMatches( xmlParser,
                                     "MaxInnerLoops",
                                     maxInnerLoops );
      InterpretElementIfNameMatches( xmlParser,
                                     "MaxOuterLoops",
                                     maxOuterLoops );
    }
    CheckSurvivalProbabilityThreshold( survivalProbabilityThreshold );
    return new CosmoTransitionsRunner(
                              InterpretTunnelingStrategy( tunnelingStrategy ),
                                       survivalProbabilityThreshold,
                                       temperatureAccuracy,
                                       pathToCosmotransitions,
                                       resolutionOfDsbVacuum,
                                       maxInnerLoops,
                                       maxOuterLoops,
                                       thermalStraightPathFitResolution,
                                       vacuumSeparationFraction );
  }

  // This parses arguments from constructorArguments and uses them to
  // construct a BounceAlongPathWithThreshold instance to use to calculate
  // the tunneling from the DSB vacuum to the panic vacuum, if possible.
  TunnelingCalculator* VevaciousPlusPlus::SetUpBounceAlongPathWithThreshold(
                                      std::string const& constructorArguments )
  {
    std::string tunnelPathFinders( "" );
    std::string bouncePotentialFitClass( "BubbleShootingOnSpline" );
    std::string bouncePotentialFitArguments( "" );
    std::string tunnelingStrategy( "ThermalThenQuantum" );
    double survivalProbabilityThreshold( 0.1 );
    int thermalIntegrationResolution( 5 );
    int temperatureAccuracy( 7 );
    int resolutionOfPathPotential( 100 );
    double vacuumSeparationFraction( 0.2 );

    BOL::AsciiXmlParser xmlParser;
    xmlParser.loadString( constructorArguments );
    while( xmlParser.readNextElement() )
    {
      InterpretElementIfNameMatches( xmlParser,
                                     "TunnelingStrategy",
                                     tunnelingStrategy );
      InterpretElementIfNameMatches( xmlParser,
                                     "SurvivalProbabilityThreshold",
                                     survivalProbabilityThreshold );
      InterpretElementIfNameMatches( xmlParser,
                                     "ThermalActionResolution",
                                     thermalIntegrationResolution );
      InterpretElementIfNameMatches( xmlParser,
                                     "CriticalTemperatureAccuracy",
                                     temperatureAccuracy );
      InterpretElementIfNameMatches( xmlParser,
                                     "PathResolution",
                                     resolutionOfPathPotential );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinimumVacuumSeparationFraction",
                                     vacuumSeparationFraction );
      ReadClassAndArguments( xmlParser,
                             "BouncePotentialFit",
                             bouncePotentialFitClass,
                             bouncePotentialFitArguments );
      InterpretElementIfNameMatches( xmlParser,
                                     "TunnelPathFinders",
                                     tunnelPathFinders );
    }
    CheckSurvivalProbabilityThreshold( survivalProbabilityThreshold );
    std::vector< BouncePathFinder* > pathFinders;
    SetUpBouncePathFinders( tunnelPathFinders,
                            pathFinders );
    if( pathFinders.empty() )
    {
      std::stringstream errorStream;
      errorStream
      << "<TunnelPathFinders> produced no valid tunnel path finding objects!";
      throw std::runtime_error( errorStream.str() );
    }
    return new BounceAlongPathWithThreshold( pathFinders,
                          SetUpBounceActionCalculator( bouncePotentialFitClass,
                                                 bouncePotentialFitArguments ),
                               InterpretTunnelingStrategy( tunnelingStrategy ),
                                             survivalProbabilityThreshold,
                                             thermalIntegrationResolution,
                                             temperatureAccuracy,
                                             resolutionOfPathPotential,
                                             vacuumSeparationFraction );
  }

  // This parses the XMl of tunnelPathFinders to construct a set of
  // BouncePathFinder instances, filling pathFinders with pointers to them.
  void VevaciousPlusPlus::SetUpBouncePathFinders(
                                          std::string const& tunnelPathFinders,
                                std::vector< BouncePathFinder* >& pathFinders )
  {
    BOL::AsciiXmlParser xmlParser;
    xmlParser.loadString( tunnelPathFinders );
    std::string classType( "" );
    std::string constructorArguments( "" );
    while( xmlParser.readNextElement() )
    {
      classType.clear();
      constructorArguments.clear();
      ReadClassAndArguments( xmlParser,
                             "PathFinder",
                             classType,
                             constructorArguments );
      if( !( classType.empty()
             ||
             constructorArguments.empty() ) )
      {
        if( classType.compare( "MinuitOnPotentialOnParallelPlanes" ) == 0 )
        {
          pathFinders.push_back( CreateMinuitOnPotentialOnParallelPlanes(
                                                      constructorArguments ) );
        }
        else if( classType.compare( "MinuitOnPotentialPerpendicularToPath" )
                 == 0 )
        {
          pathFinders.push_back( CreateMinuitOnPotentialPerpendicularToPath(
                                                      constructorArguments ) );
        }
        else if( classType.compare( "MinuitOnPathPerpendicularForces" )
                 == 0 )
        {
          pathFinders.push_back( CreateMinuitOnPathNormalInertialPotential(
                                                      constructorArguments ) );
        }
      }
    }
  }

  // This parses arguments from constructorArguments and uses them to
  // construct a MinuitOnPotentialOnParallelPlanes instance to use to try to
  // extremize the bounce action.
  BouncePathFinder* VevaciousPlusPlus::CreateMinuitOnPotentialOnParallelPlanes(
                                      std::string const& constructorArguments )
  {
    // The <ConstructorArguments> for this class should have child elements
    // <NumberOfPathSegments>, <MinuitStrategy> and <MinuitTolerance>.
    unsigned int numberOfPathSegments( 100 );
    unsigned int minuitStrategy( 1 );
    double minuitToleranceFraction( 0.5 );
    BOL::AsciiXmlParser xmlParser;
    xmlParser.loadString( constructorArguments );
    while( xmlParser.readNextElement() )
    {
      InterpretElementIfNameMatches( xmlParser,
                                     "NumberOfPathSegments",
                                     numberOfPathSegments );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinuitStrategy",
                                     minuitStrategy );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinuitTolerance",
                                     minuitToleranceFraction );
    }
    return new MinuitOnPotentialOnParallelPlanes( numberOfPathSegments,
                                                  minuitStrategy,
                                                  minuitToleranceFraction );
  }

  // This parses arguments from constructorArguments and uses them to
  // construct a MinuitOnPotentialPerpendicularToPath instance to use to try to
  // extremize the bounce action.
  BouncePathFinder*
  VevaciousPlusPlus::CreateMinuitOnPotentialPerpendicularToPath(
                                      std::string const& constructorArguments )
  {
    // The <ConstructorArguments> for this class should have child elements
    // <NumberOfPathSegments>, <MinuitStrategy> and <MinuitTolerance>.
    int numberOfPathSegments( 100 );
    int numberOfAllowedWorsenings( 3 );
    double convergenceThresholdFraction( 0.05 );
    double minuitDampingFraction( 0.75 );
    std::string neighborDisplacementWeightsString( "0.5, 0.25" );
    int minuitStrategy( 1 );
    double minuitToleranceFraction( 0.5 );
    BOL::AsciiXmlParser xmlParser;
    xmlParser.loadString( constructorArguments );
    while( xmlParser.readNextElement() )
    {
      InterpretElementIfNameMatches( xmlParser,
                                     "NumberOfPathSegments",
                                     numberOfPathSegments );
      InterpretElementIfNameMatches( xmlParser,
                                     "NumberOfAllowedWorsenings",
                                     numberOfAllowedWorsenings );
      InterpretElementIfNameMatches( xmlParser,
                                     "ConvergenceThresholdFraction",
                                     convergenceThresholdFraction );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinuitDampingFraction",
                                     minuitDampingFraction );
      InterpretElementIfNameMatches( xmlParser,
                                     "NeighborDisplacementWeights",
                                     neighborDisplacementWeightsString );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinuitStrategy",
                                     minuitStrategy );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinuitTolerance",
                                     minuitToleranceFraction );
    }
    BOL::VectorlikeArray< std::string > weightsStrings;
    BOL::StringParser::parseByChar( neighborDisplacementWeightsString,
                                    weightsStrings,
                                    " \t\n,;" );
    std::vector< double > neighborDisplacementWeights;
    for( int weightIndex( 0 );
         weightIndex < weightsStrings.getSize();
         ++weightIndex )
    {
      neighborDisplacementWeights.push_back(
          BOL::StringParser::stringToDouble( weightsStrings[ weightIndex ] ) );
    }

    return new MinuitOnPotentialPerpendicularToPath( numberOfPathSegments,
                                                     numberOfAllowedWorsenings,
                                                  convergenceThresholdFraction,
                                                     minuitDampingFraction,
                                                   neighborDisplacementWeights,
                                                     minuitStrategy,
                                                     minuitToleranceFraction );
  }

  // This parses arguments from constructorArguments and uses them to
  // construct a MinuitOnPathNormalInertialPotential instance to use to try to
  // extremize the bounce action.
  BouncePathFinder*
  VevaciousPlusPlus::CreateMinuitOnPathNormalInertialPotential(
                                      std::string const& constructorArguments )
  {
    // The <ConstructorArguments> for this class should have child elements
    // <NumberOfPathSegments>, <MinuitStrategy> and <MinuitTolerance>.
    int numberOfPathSegments( 100 );
    int numberOfAllowedWorsenings( 3 );
    double convergenceThresholdFraction( 0.05 );
    double minuitDampingFraction( 0.75 );
    std::string neighborDisplacementWeightsString( "0.5, 0.25" );
    int minuitStrategy( 1 );
    double minuitToleranceFraction( 0.5 );
    BOL::AsciiXmlParser xmlParser;
    xmlParser.loadString( constructorArguments );
    while( xmlParser.readNextElement() )
    {
      InterpretElementIfNameMatches( xmlParser,
                                     "NumberOfPathSegments",
                                     numberOfPathSegments );
      InterpretElementIfNameMatches( xmlParser,
                                     "NumberOfAllowedWorsenings",
                                     numberOfAllowedWorsenings );
      InterpretElementIfNameMatches( xmlParser,
                                     "ConvergenceThresholdFraction",
                                     convergenceThresholdFraction );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinuitDampingFraction",
                                     minuitDampingFraction );
      InterpretElementIfNameMatches( xmlParser,
                                     "NeighborDisplacementWeights",
                                     neighborDisplacementWeightsString );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinuitStrategy",
                                     minuitStrategy );
      InterpretElementIfNameMatches( xmlParser,
                                     "MinuitTolerance",
                                     minuitToleranceFraction );
    }
    BOL::VectorlikeArray< std::string > weightsStrings;
    BOL::StringParser::parseByChar( neighborDisplacementWeightsString,
                                    weightsStrings,
                                    " \t\n,;" );
    std::vector< double > neighborDisplacementWeights;
    for( int weightIndex( 0 );
         weightIndex < weightsStrings.getSize();
         ++weightIndex )
    {
      neighborDisplacementWeights.push_back(
          BOL::StringParser::stringToDouble( weightsStrings[ weightIndex ] ) );
    }

    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "CreateMinuitOnPathNormalInertialPotential(...) returning NULL as the"
    << " class does not yet exist!";
    std::cout << std::endl;
    return NULL;/**/
    /*return new MinuitOnPathNormalInertialPotential( numberOfPathSegments,
                                                     numberOfAllowedWorsenings,
                                                  convergenceThresholdFraction,
                                                     minuitDampingFraction,
                                                   neighborDisplacementWeights,
                                                     minuitStrategy,
                                                   minuitToleranceFraction );*/
  }

  // This parses arguments from constructorArguments and uses them to
  // construct a BubbleShootingOnPathInFieldSpace instance to use to calculate
  // the bounce action on given paths.
  BounceActionCalculator*
  VevaciousPlusPlus::SetUpBubbleShootingOnPathInFieldSpace(
                                      std::string const& constructorArguments )
  {
    double lengthScaleResolutionForBounce( 0.05 );
    unsigned int shootAttemptsForBounce( 32 );

    BOL::AsciiXmlParser xmlParser;
    xmlParser.loadString( constructorArguments );
    while( xmlParser.readNextElement() )
    {
      InterpretElementIfNameMatches( xmlParser,
                                     "RadialResolution",
                                     lengthScaleResolutionForBounce );
      InterpretElementIfNameMatches( xmlParser,
                                     "NumberShootAttemptsAllowed",
                                     shootAttemptsForBounce );
    }

    return
    new BubbleShootingOnPathInFieldSpace( lengthScaleResolutionForBounce,
                                          shootAttemptsForBounce );
  }

} /* namespace VevaciousPlusPlus */
