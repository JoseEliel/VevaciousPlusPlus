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
  // in that file. Since this constructor leaves deleterForPotentialFunction,
  // deleterForPotentialMinimizer, and deleterForPotentialMinimizer all as
  // NULL, and no other function sets them, there should be no problem with
  // the destructor calling delete on these pointers, as they do not get set to
  // point at the addresses of the given components.
  VevaciousPlusPlus::VevaciousPlusPlus( SlhaManager& slhaManager,
                                        PotentialMinimizer& potentialMinimizer,
                                   TunnelingCalculator& tunnelingCalculator ) :
    potentialFunction( NULL ),
    deleterForPotentialFunction( NULL ),
    runningParameterManager(),
    potentialMinimizer( &potentialMinimizer ),
    deleterForPotentialMinimizer( NULL ),
    tunnelingCalculator( &tunnelingCalculator ),
    deleterForTunnelingCalculator( NULL ),
    currentTime()
  {
    // This constructor is just an initialization list.
  }

  VevaciousPlusPlus::VevaciousPlusPlus(
                                  std::string const& initializationFileName ) :
    potentialFunction( NULL ),
    deleterForPotentialFunction( NULL ),
    runningParameterManager(),
    potentialMinimizer( NULL ),
    deleterForPotentialMinimizer( NULL ),
    tunnelingCalculator( NULL ),
    deleterForTunnelingCalculator( NULL ),
    currentTime()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "NEED TO ENSURE THAT ALL HEADERS ARE PROPERLY INCLUDED!";
    std::cout << std::endl;/**/
    PotentialFromPolynomialAndMasses* potentialFromPolynomialAndMasses( NULL );
    std::string potentialClass( "FixedScaleOneLoopPotential" );
    std::string potentialArguments( "" );
    std::string minimizerClass( "GradientFromStartingPoints" );
    std::string minimizerArguments( "" );
    std::string tunnelingClass( "BounceAlongPathWithThreshold" );
    std::string tunnelingArguments( "" );
    BOL::AsciiXmlParser fileParser;
    fileParser.openRootElementOfFile( initializationFileName );
    while( fileParser.readNextElement() )
    {
      if( fileParser.currentElementNameMatches( "PotentialClass" ) )
      {
        BOL::AsciiXmlParser oneNestedParser;
        oneNestedParser.loadString(
                                fileParser.getTrimmedCurrentElementContent() );
        while( oneNestedParser.readNextElement() )
        {
          if( oneNestedParser.currentElementNameMatches( "ClassType" ) )
          {
            potentialClass.assign(
                           oneNestedParser.getTrimmedCurrentElementContent() );
          }
          else if( oneNestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
          {
            potentialArguments.assign(
                         oneNestedParser.getTrimmedCurrentElementContent() );
          }
        }
      }
      else if( fileParser.currentElementNameMatches( "MinimizerClass" ) )
      {
        BOL::AsciiXmlParser oneNestedParser;
        oneNestedParser.loadString(
                                fileParser.getTrimmedCurrentElementContent() );
        while( oneNestedParser.readNextElement() )
        {
          if( oneNestedParser.currentElementNameMatches( "ClassType" ) )
          {
            minimizerClass.assign(
                           oneNestedParser.getTrimmedCurrentElementContent() );
          }
          else if( oneNestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
          {
            minimizerArguments.assign(
                         oneNestedParser.getTrimmedCurrentElementContent() );
          }
        }
      }
      else if( fileParser.currentElementNameMatches( "TunnelingClass" ) )
      {
        BOL::AsciiXmlParser oneNestedParser;
        oneNestedParser.loadString(
                                fileParser.getTrimmedCurrentElementContent() );
        while( oneNestedParser.readNextElement() )
        {
          if( oneNestedParser.currentElementNameMatches( "ClassType" ) )
          {
            tunnelingClass.assign(
                           oneNestedParser.getTrimmedCurrentElementContent() );
          }
          else if( oneNestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
          {
            tunnelingArguments.assign(
                         oneNestedParser.getTrimmedCurrentElementContent() );
          }
        }
      }
    }

    // Now we know the overall picture of what components to set up.
    // ALL SUB-COMPONENTS FOR potentialMinimizer AND tunnelingCalculator ARE
    // MEMORY-MANAGED BY THE COMPONENTS THEMSELVES! The VevaciousPlusPlus
    // destructor only deletes deleterForTunnelingCalculator,
    // deleterForPotentialMinimizer, and deleterForPotentialFunction, while
    // this constructor is allocating much more memory than that.

    // First the potential function:
    if( potentialClass.compare( "FixedScaleOneLoopPotential" ) == 0 )
    {
      potentialFromPolynomialAndMasses
      = new FixedScaleOneLoopPotential( potentialArguments,
                                        runningParameterManager );
    }
    else if( potentialClass.compare( "RgeImprovedOneLoopPotential" ) == 0 )
    {
      potentialFromPolynomialAndMasses
      = new RgeImprovedOneLoopPotential( potentialArguments,
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
    potentialFunction = (PotentialFunction*)potentialFromPolynomialAndMasses;
    deleterForPotentialFunction = potentialFunction;

    // Next the potential minimizer:
    if( minimizerClass.compare( "GradientFromStartingPoints" ) == 0 )
    {
      // We need to assemble the components for a GradientFromStartingPoints
      // object: a StartingPointFinder and a GradientMinimizer. We need to find
      // out what derived classes to actually use.
      BOL::AsciiXmlParser elementParser;
      elementParser.loadString( minimizerArguments );
      StartingPointFinder* startingPointFinder( NULL );
      std::string startingPointFinderClass( "Hom4ps2Runner" );
      std::string startingPointFinderArguments( "" );
      GradientMinimizer* gradientMinimizer( NULL );
      std::string gradientMinimizerClass( "MinuitPotentialMinimizer" );
      std::string gradientMinimizerArguments( "" );
      double extremumSeparationThresholdFraction( 0.05 );
      double nonDsbRollingToDsbScalingFactor( 10.0 );

      // The <ConstructorArguments> for this class should have child elements
      // <StartingPointFinderClass> and <GradientMinimizerClass>, and
      // optionally <ExtremumSeparationThresholdFraction> and
      // <NonDsbRollingToDsbScalingFactor>.
      while( elementParser.readNextElement() )
      {
        if( elementParser.currentElementNameMatches(
                                                 "StartingPointFinderClass" ) )
        {
          // <StartingPointFinderClass> should have child elements <ClassName>
          // and <ConstructorArguments>.
          BOL::AsciiXmlParser nestedParser;
          nestedParser.loadString(
                             elementParser.getTrimmedCurrentElementContent() );
          while( nestedParser.readNextElement() )
          {
            if( nestedParser.currentElementNameMatches( "ClassType" ) )
            {
              startingPointFinderClass.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
            }
            else if( nestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
            {
              startingPointFinderArguments.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
            }
          }
        }
        else if( elementParser.currentElementNameMatches(
                                                   "GradientMinimizerClass" ) )
        {
          // <GradientMinimizerClass> should have child elements <ClassName>
          // and <ConstructorArguments>.
          BOL::AsciiXmlParser nestedParser;
          nestedParser.loadString(
                             elementParser.getTrimmedCurrentElementContent() );
          while( nestedParser.readNextElement() )
          {
            if( nestedParser.currentElementNameMatches( "ClassType" ) )
            {
              gradientMinimizerClass.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
            }
            else if( nestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
            {
              gradientMinimizerArguments.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
            }
          }
        }
        else if( elementParser.currentElementNameMatches(
                                      "ExtremumSeparationThresholdFraction" ) )
        {
          extremumSeparationThresholdFraction
          = BOL::StringParser::stringToDouble(
                           elementParser.getTrimmedCurrentElementContent() );
        }
        else if( elementParser.currentElementNameMatches(
                                          "NonDsbRollingToDsbScalingFactor" ) )
        {
          nonDsbRollingToDsbScalingFactor
          = BOL::StringParser::stringToDouble(
                           elementParser.getTrimmedCurrentElementContent() );
        }
      }

      if( startingPointFinderClass.compare( "Hom4ps2Runner" ) == 0 )
      {
        startingPointFinder = new Hom4ps2Runner(
          potentialFromPolynomialAndMasses->HomotopyContinuationTargetSystem(),
                                              startingPointFinderArguments );
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
        gradientMinimizer = new MinuitPotentialMinimizer( *potentialFunction,
                                                  gradientMinimizerArguments );
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
      potentialMinimizer
      = new GradientFromStartingPoints( startingPointFinder,
                                        gradientMinimizer,
                                        extremumSeparationThresholdFraction,
                                        nonDsbRollingToDsbScalingFactor );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<MinimizerClass> was not a recognized form! The only type"
      << " currently valid is \"GradientFromStartingPoints\".";
      throw std::runtime_error( errorStream.str() );
    }
    deleterForPotentialMinimizer = potentialMinimizer;

    // Next the tunneling calculator:
    if( tunnelingClass.compare( "CosmoTransitionsRunner" ) == 0 )
    {
      tunnelingCalculator
      = new CosmoTransitionsRunner( *potentialFromPolynomialAndMasses,
                                    *potentialFunction,
                                    tunnelingArguments );
    }
    else if( tunnelingClass.compare( "BounceAlongPathWithThreshold" ) == 0 )
    {
      // We need to assemble the components for a BounceAlongPathWithThreshold
      // object: a BounceActionCalculator and a BouncePathFinder. We need to
      // find out what derived classes to actually use.
      BOL::AsciiXmlParser elementParser;
      elementParser.loadString( tunnelingArguments );
      BounceActionCalculator* bounceActionCalculator( NULL );
      std::string bouncePotentialFitClass( "BubbleShootingOnSpline" );
      std::string bouncePotentialFitArguments( "" );
      BouncePathFinder* bouncePathFinder( NULL );
      std::string tunnelPathFinderClass( "MinuitNodePotentialMinimizer" );
      std::string tunnelPathFinderArguments( "" );

      // The <ConstructorArguments> for this class should have child elements
      // <BouncePotentialFit> and <TunnelPathFinder>.
      while( elementParser.readNextElement() )
      {
        if( elementParser.currentElementNameMatches( "BouncePotentialFit" ) )
        {
          // <BouncePotentialFit> should have child elements <ClassName> and
          // <ConstructorArguments>.
          BOL::AsciiXmlParser nestedParser;
          nestedParser.loadString(
                             elementParser.getTrimmedCurrentElementContent() );
          while( nestedParser.readNextElement() )
          {
            if( nestedParser.currentElementNameMatches( "ClassType" ) )
            {
              bouncePotentialFitClass.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
            }
            else if( nestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
            {
              bouncePotentialFitArguments.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
            }
          }
        }
        else if( elementParser.currentElementNameMatches(
                                                         "TunnelPathFinder" ) )
        {
          // <TunnelPathFinder> should have child elements <ClassName> and
          // <ConstructorArguments>.
          BOL::AsciiXmlParser nestedParser;
          nestedParser.loadString(
                           elementParser.getTrimmedCurrentElementContent() );
          while( nestedParser.readNextElement() )
          {
            if( nestedParser.currentElementNameMatches( "ClassType" ) )
            {
              tunnelPathFinderClass.assign(
                           nestedParser.getTrimmedCurrentElementContent() );
            }
            else if( nestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
            {
              tunnelPathFinderArguments.assign(
                           nestedParser.getTrimmedCurrentElementContent() );
            }
          }
        }
      }

      if( bouncePotentialFitClass.compare( "BubbleShootingOnSpline" ) == 0 )
      {
        bounceActionCalculator = new BubbleShootingOnSpline( potentialFunction,
                                                 bouncePotentialFitArguments );
      }
      else
      {
        std::stringstream errorStream;
        errorStream
        << "<BouncePotentialFit> was not a recognized form! The only type"
        << " currently valid is \"BubbleShootingOnSpline\".";
        throw std::runtime_error( errorStream.str() );
      }

      if( ( tunnelPathFinderClass.compare( "MinuitNodePotentialMinimizer" )
            == 0 )
          ||
          ( tunnelPathFinderClass.compare( "MinuitPathBounceMinimizer" )
            == 0 )
          ||
          ( tunnelPathFinderClass.compare( "MinuitPathPotentialMinimizer" )
            == 0 ) )
      {
        size_t numberOfPotentialSamplePoints( 15 );
        size_t movesPerImprovement( 100 );
        unsigned int minuitStrategy( 1 );
        double minuitTolerance( 0.5 );
        std::string pathParameterizationClass( "PathFromNodes" );
        std::string pathParameterizationArguments( "" );
        std::string pathBouncePotentialFitClass( "BubbleShootingOnSpline" );
        std::string pathBouncePotentialFitArguments( "" );

        BOL::AsciiXmlParser nestedParser;
        nestedParser.loadString(
                             elementParser.getTrimmedCurrentElementContent() );
        while( nestedParser.readNextElement() )
        {
          if( nestedParser.currentElementNameMatches(
                                                      "MovesPerImprovement" ) )
          {
            movesPerImprovement = BOL::StringParser::stringToInt(
                              nestedParser.getTrimmedCurrentElementContent() );
          }
          else if( nestedParser.currentElementNameMatches( "MinuitStrategy" ) )
          {
            minuitStrategy = BOL::StringParser::stringToInt(
                              nestedParser.getTrimmedCurrentElementContent() );
          }
          else if( nestedParser.currentElementNameMatches(
                                                          "MinuitTolerance" ) )
          {
            minuitTolerance = BOL::StringParser::stringToDouble(
                              nestedParser.getTrimmedCurrentElementContent() );
          }
          else if( nestedParser.currentElementNameMatches(
                                            "NumberOfPotentialSamplePoints" ) )
          {
            numberOfPotentialSamplePoints = BOL::StringParser::stringToInt(
                              nestedParser.getTrimmedCurrentElementContent() );
          }
          else if( nestedParser.currentElementNameMatches(
                                                       "BouncePotentialFit" ) )
          {
            BOL::AsciiXmlParser doublyNestedParser;
            doublyNestedParser.loadString(
                              nestedParser.getTrimmedCurrentElementContent() );
            while( doublyNestedParser.readNextElement() )
            {
              if( doublyNestedParser.currentElementNameMatches( "ClassType" ) )
              {
                pathBouncePotentialFitClass.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
              }
              else if( doublyNestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
              {
                pathBouncePotentialFitArguments.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
              }
            }
          }
          else if( nestedParser.currentElementNameMatches(
                                                     "PathParameterization" ) )
          {
            BOL::AsciiXmlParser doublyNestedParser;
            doublyNestedParser.loadString(
                              nestedParser.getTrimmedCurrentElementContent() );
            while( doublyNestedParser.readNextElement() )
            {
              if( doublyNestedParser.currentElementNameMatches( "ClassType" ) )
              {
                pathParameterizationClass.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
              }
              else if( doublyNestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
              {
                pathParameterizationArguments.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
              }
            }
          }
        }

        // We still need to extract the type of path parameterization.
        nestedParser.loadString( pathParameterizationArguments );
        if( pathParameterizationClass.compare(
                                      "PathFromPolynomialCoefficients" ) == 0 )
        {
          // placeholder:
          /**/std::cout << std::endl
          << "Placeholder: "
          << "DO SOMETHING HERE! NULL POINTER ABOUT TO CAUSE SEG FAULT!";
          std::cout << std::endl;/**/
        }
        bouncePathFinder = NULL;
        // bouncePathFinder = new MinuitNodePotentialMinimizer( potentialFunction,
        //                                           tunnelPathFinderArguments );
      }
      else
      {
        std::stringstream errorStream;
        errorStream
        << "<TunnelPathFinder> was not a recognized form! The only types"
        << " currently valid are \"MinuitNodePotentialMinimizer\","
        << " \"MinuitPathBounceMinimizer\", or"
        << " \"MinuitPathPotentialMinimizer\".";
        throw std::runtime_error( errorStream.str() );
      }

      // Now we have the components for tunnelingCalculator:
      tunnelingCalculator
      = new BounceAlongPathWithThreshold( (*potentialFunction),
                                          bounceActionCalculator,
                                          bouncePathFinder,
                                          tunnelingArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<TunnelingClass> was not a recognized form! The only types"
      << " currently valid are \"BounceAlongPathWithThreshold\" or"
      << " \"CosmoTransitionsRunner\".";
      throw std::runtime_error( errorStream.str() );
    }
    deleterForTunnelingCalculator = tunnelingCalculator;
  }

  VevaciousPlusPlus::~VevaciousPlusPlus()
  {
    delete deleterForTunnelingCalculator;
    delete deleterForPotentialMinimizer;
    delete deleterForPotentialFunction;
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
