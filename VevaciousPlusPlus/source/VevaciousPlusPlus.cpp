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
    slhaManager( &slhaManager ),
    deleterForSlhaManager( NULL ),
    potentialFunction( NULL ),
    deleterForPotentialFunction( NULL ),
    potentialMinimizer( &potentialMinimizer ),
    deleterForPotentialMinimizer( NULL ),
    tunnelingCalculator( &tunnelingCalculator ),
    deleterForTunnelingCalculator( NULL ),
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
    slhaManager( NULL ),
    deleterForSlhaManager( NULL ),
    potentialFunction( NULL ),
    deleterForPotentialFunction( NULL ),
    potentialMinimizer( NULL ),
    deleterForPotentialMinimizer( NULL ),
    tunnelingCalculator( NULL ),
    deleterForTunnelingCalculator( NULL ),
    currentTime()
  {
    RunningParameterManager*
    runningParameterManager( new RunningParameterManager() );
    slhaManager = runningParameterManager;
    deleterForSlhaManager = slhaManager;
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
        BOL::AsciiXmlParser elementParser;
        elementParser.loadString(
                                fileParser.getTrimmedCurrentElementContent() );
        while( elementParser.readNextElement() )
        {
          if( elementParser.currentElementNameMatches( "ClassType" ) )
          {
            potentialClass.assign(
                             elementParser.getTrimmedCurrentElementContent() );
          }
          else if( elementParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
          {
            potentialArguments.assign(
                             elementParser.getTrimmedCurrentElementContent() );
          }
        }
      }
      else if( fileParser.currentElementNameMatches( "MinimizerClass" ) )
      {
        BOL::AsciiXmlParser elementParser;
        elementParser.loadString(
                                fileParser.getTrimmedCurrentElementContent() );
        while( elementParser.readNextElement() )
        {
          if( elementParser.currentElementNameMatches( "ClassType" ) )
          {
            minimizerClass.assign(
                             elementParser.getTrimmedCurrentElementContent() );
          }
          else if( elementParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
          {
            minimizerArguments.assign(
                             elementParser.getTrimmedCurrentElementContent() );
          }
        }
      }
      else if( fileParser.currentElementNameMatches( "TunnelingClass" ) )
      {
        BOL::AsciiXmlParser elementParser;
        elementParser.loadString(
                                fileParser.getTrimmedCurrentElementContent() );
        while( elementParser.readNextElement() )
        {
          if( elementParser.currentElementNameMatches( "ClassType" ) )
          {
            tunnelingClass.assign(
                             elementParser.getTrimmedCurrentElementContent() );
          }
          else if( elementParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
          {
            tunnelingArguments.assign(
                             elementParser.getTrimmedCurrentElementContent() );
          }
        }
      }
    }

    // Now we know the overall picture of what components to set up.
    // ALL SUB-COMPONENTS FOR potentialMinimizer AND tunnelingCalculator ARE
    // MEMORY-MANAGED BY THE COMPONENTS THEMSELVES! The VevaciousPlusPlus
    // destructor only deletes deleterForTunnelingCalculator,
    // deleterForPotentialMinimizer, and deleterForPotentialFunction, while
    // this constructor is allocating much more memory than that. Everything
    // would be clear if we restricted ourselves to requiring a C++11-compliant
    // compiler, as then we could have the constructors use std::unique_ptrs,
    // but alas, we're sticking to C++98.

    // First the potential function: since the only options are both derived
    // from PotentialFromPolynomialAndMasses with no specific extra arguments,
    // we extract them from the XML here.


    BOL::AsciiXmlParser elementParser;
    elementParser.loadString( potentialArguments );
    std::string modelFilename( "./ModelFiles/SM.vin" );
    double scaleRangeMinimumFactor( 10.0 );
    bool treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions( false );
    while( elementParser.readNextElement() )
    {
      if( elementParser.currentElementNameMatches( "ModelFile" ) )
      {
        modelFilename.assign(
                             elementParser.getTrimmedCurrentElementContent() );
      }
      else if( elementParser.currentElementNameMatches( "RollOnlyMinima" ) )
      {
        std::string rollOnlyMinima(
                             elementParser.getTrimmedCurrentElementContent() );
        BOL::StringParser::transformToLowercase( rollOnlyMinima );
        treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions
        = ( rollOnlyMinima.compare( "true" ) == 0 );
        if( !treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions
            &&
            ( rollOnlyMinima.compare( "yes" ) == 0 ) )
        {
          treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions = true;
        }
      }
      else if( elementParser.currentElementNameMatches(
                                                  "ScaleRangeMinimumFactor" ) )
      {
        scaleRangeMinimumFactor = BOL::StringParser::stringToDouble(
                             elementParser.getTrimmedCurrentElementContent() );
      }
    }

    if( potentialClass.compare( "FixedScaleOneLoopPotential" ) == 0 )
    {
      potentialFromPolynomialAndMasses
      = new FixedScaleOneLoopPotential( modelFilename,
                                        scaleRangeMinimumFactor,
                       treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions,
                                        *runningParameterManager );
    }
    else if( potentialClass.compare( "RgeImprovedOneLoopPotential" ) == 0 )
    {
      potentialFromPolynomialAndMasses
      = new RgeImprovedOneLoopPotential( modelFilename,
                                         scaleRangeMinimumFactor,
                       treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions,
                                         *runningParameterManager );
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
    potentialFunction = potentialFromPolynomialAndMasses;
    deleterForPotentialFunction = potentialFunction;

    size_t const numberOfFields( potentialFunction->NumberOfFieldVariables() );

    // Next the potential minimizer:
    if( minimizerClass.compare( "GradientFromStartingPoints" ) == 0 )
    {
      // We need to assemble the components for a GradientFromStartingPoints
      // object: a StartingPointFinder and a GradientMinimizer. We need to find
      // out what derived classes to actually use.
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
        elementParser.loadString( startingPointFinderArguments );
        // The <ConstructorArguments> for this class should have child elements
        // <PathToHom4ps2> and <Hom4ps2Argument>.
        std::string pathToHom4ps2( "./HOM4PS2/" );
        std::string homotopyType( "2" );
        while( elementParser.readNextElement() )
        {
          if( elementParser.currentElementNameMatches( "PathToHom4ps2" ) )
          {
            pathToHom4ps2.assign(
                             elementParser.getTrimmedCurrentElementContent() );
          }
          else if( elementParser.currentElementNameMatches(
                                                          "Hom4ps2Argument" ) )
          {
            homotopyType.assign(
                             elementParser.getTrimmedCurrentElementContent() );
          }
        }
        startingPointFinder = new Hom4ps2Runner(
                      *(potentialFunction->HomotopyContinuationTargetSystem()),
                                                 pathToHom4ps2,
                                                 homotopyType );
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
        elementParser.loadString( gradientMinimizerArguments );
        double errorFraction( 0.1 );
        double errorMinimum( 1.0 );
        unsigned int minuitStrategy( 1 );
        while( elementParser.readNextElement() )
        {
          if( elementParser.currentElementNameMatches(
                                                  "InitialStepSizeFraction" ) )
          {
            errorFraction = BOL::StringParser::stringToDouble(
                            elementParser.getTrimmedCurrentElementContent() );
          }
          else if( elementParser.currentElementNameMatches(
                                                   "MinimumInitialStepSize" ) )
          {
            errorMinimum = BOL::StringParser::stringToDouble(
                            elementParser.getTrimmedCurrentElementContent() );
          }
          else if( elementParser.currentElementNameMatches(
                                                           "MinuitStrategy" ) )
          {
            minuitStrategy = BOL::StringParser::stringToInt(
                            elementParser.getTrimmedCurrentElementContent() );
          }
        }
        gradientMinimizer = new MinuitPotentialMinimizer( *potentialFunction,
                                                          errorFraction,
                                                          errorMinimum,
                                                          minuitStrategy );
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
      potentialMinimizer = new GradientFromStartingPoints( *potentialFunction,
                                                           startingPointFinder,
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
    if( ( tunnelingClass.compare( "CosmoTransitionsRunner" ) == 0 )
        ||
        ( tunnelingClass.compare( "BounceAlongPathWithThreshold" ) == 0 ) )
    {
      elementParser.loadString( tunnelingArguments );
      TunnelingCalculator::TunnelingStrategy
      tunnelingStrategy( TunnelingCalculator::JustThermal );
      double survivalProbabilityThreshold( 0.1 );
      size_t temperatureAccuracy( 7 );
      size_t evaporationResolution( 3 );
      std::string pathToCosmotransitions( "./cosmoTransitions/" );
      size_t resolutionOfDsbVacuum( 20 );
      size_t maxInnerLoops( 10 );
      size_t maxOuterLoops( 10 );
      std::string bouncePotentialFitClass( "BubbleShootingOnSpline" );
      std::string bouncePotentialFitArguments( "" );
      std::string tunnelPathFinderClass( "MinuitNodePotentialMinimizer" );
      std::string tunnelPathFinderArguments( "" );
      size_t thermalIntegrationResolution( 5 );
      while( elementParser.readNextElement() )
      {
        if( elementParser.currentElementNameMatches( "TunnelingStrategy" ) )
        {
          std::string tunnelingStrategyString(
                             elementParser.getTrimmedCurrentElementContent() );
          if( ( tunnelingStrategyString.compare( "DefaultTunneling" ) == 0 )
              ||
              ( tunnelingStrategyString.compare( "ThermalThenQuantum" )
                == 0 ) )
          {
            tunnelingStrategy = TunnelingCalculator::ThermalThenQuantum;
          }
          else if( tunnelingStrategyString.compare( "QuantumThenThermal" )
                   == 0 )
          {
            tunnelingStrategy
            = TunnelingCalculator::QuantumThenThermal;
          }
          else if( tunnelingStrategyString.compare( "JustThermal" ) == 0 )
          {
            tunnelingStrategy = TunnelingCalculator::JustThermal;
          }
          else if( tunnelingStrategyString.compare( "JustQuantum" ) == 0 )
          {
            tunnelingStrategy = TunnelingCalculator::JustQuantum;
          }
          else if( ( tunnelingStrategyString.compare( "NoTunneling" ) == 0 )
                   ||
                   ( tunnelingStrategyString.compare( "None" ) == 0 ) )
          {
            tunnelingStrategy = TunnelingCalculator::NoTunneling;
          }
        }
        else if( elementParser.currentElementNameMatches(
                                             "SurvivalProbabilityThreshold" ) )
        {
          survivalProbabilityThreshold = BOL::StringParser::stringToDouble(
                             elementParser.getTrimmedCurrentElementContent() );
          if( !( ( survivalProbabilityThreshold > 0.0 )
                 &&
                 ( survivalProbabilityThreshold < 1.0 ) ) )
          {
            std::stringstream errorStream;
            errorStream
            << "<SurvivalProbabilityThreshold> must be greater than 0.0 and"
            << " less than 1.0 to be valid! Instead, the XML element gave \""
            << elementParser.getTrimmedCurrentElementContent() << "\".";
            throw std::runtime_error( errorStream.str() );
          }
        }
        else if( elementParser.currentElementNameMatches(
                                              "CriticalTemperatureAccuracy" ) )
        {
          temperatureAccuracy = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
        }
        else if( elementParser.currentElementNameMatches(
                                             "EvaporationBarrierResolution" ) )
        {
          evaporationResolution = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
        }
        else if( elementParser.currentElementNameMatches(
                                                   "PathToCosmotransitions" ) )
        {
          pathToCosmotransitions.assign(
                             elementParser.getTrimmedCurrentElementContent() );
        }
        else if( elementParser.currentElementNameMatches( "PathResolution" ) )
        {
          resolutionOfDsbVacuum = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
        }
        else if( elementParser.currentElementNameMatches( "MaxInnerLoops" ) )
        {
          maxInnerLoops = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
        }
        else if( elementParser.currentElementNameMatches( "MaxOuterLoops" ) )
        {
          maxOuterLoops = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
        }
        else if( elementParser.currentElementNameMatches(
                                                       "BouncePotentialFit" ) )
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
        else if( elementParser.currentElementNameMatches(
                                             "ThermalIntegrationResolution" ) )
        {
          thermalIntegrationResolution
          = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
        }
      }

      if( tunnelingClass.compare( "CosmoTransitionsRunner" ) == 0 )
      {
        tunnelingCalculator
        = new CosmoTransitionsRunner( *potentialFromPolynomialAndMasses,
                                      *potentialFunction,
                                      tunnelingStrategy,
                                      survivalProbabilityThreshold,
                                      temperatureAccuracy,
                                      evaporationResolution,
                                      pathToCosmotransitions,
                                      resolutionOfDsbVacuum,
                                      maxInnerLoops,
                                      maxOuterLoops );
      }
      else
      // if( tunnelingClass.compare( "BounceAlongPathWithThreshold" ) == 0 )
      // This is already true because of the outer if statement.
      {
        BounceActionCalculator* bounceActionCalculator( NULL );
        BouncePathFinder* bouncePathFinder( NULL );
        size_t numberOfPotentialSegmentsForBounce( 16 );
        size_t shootAttemptsForBounce( 32 );
        double lengthScaleResolutionForBounce( 0.05 );
        // We need to assemble the components for a
        // BounceAlongPathWithThreshold object: a BounceActionCalculator and a
        // BouncePathFinder. We need to find out what derived classes to
        // actually use.
        if( bouncePotentialFitClass.compare( "BubbleShootingOnSpline" ) == 0 )
        {
          elementParser.loadString( bouncePotentialFitArguments );
          while( elementParser.readNextElement() )
          {
            if( elementParser.currentElementNameMatches(
                                             "NumberOfNodesForPotentialFit" ) )
            {
              numberOfPotentialSegmentsForBounce
              = ( 1 + BOL::StringParser::stringToInt(
                           elementParser.getTrimmedCurrentElementContent() ) );
            }
            else if( elementParser.currentElementNameMatches(
                                               "NumberShootAttemptsAllowed" ) )
            {
              shootAttemptsForBounce = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
            }
            else if( elementParser.currentElementNameMatches(
                                                         "RadialResolution" ) )
            {
              lengthScaleResolutionForBounce
              = BOL::StringParser::stringToDouble(
                             elementParser.getTrimmedCurrentElementContent() );
            }
          }
          bounceActionCalculator
          = new BubbleShootingOnSpline( *potentialFunction,
                                        numberOfPotentialSegmentsForBounce,
                                        lengthScaleResolutionForBounce,
                                        shootAttemptsForBounce );
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
          double minuitToleranceFraction( 0.5 );
          std::string pathParameterizationClass( "PathFromNodes" );
          std::string pathParameterizationArguments( "" );
          std::string pathBouncePotentialFitClass( "BubbleShootingOnSpline" );
          std::string pathBouncePotentialFitArguments( "" );
          double nodeMovementThreshold( 0.01 );

          elementParser.loadString( tunnelPathFinderArguments );
          while( elementParser.readNextElement() )
          {
            if( elementParser.currentElementNameMatches(
                                                      "MovesPerImprovement" ) )
            {
              movesPerImprovement = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
            }
            else if( elementParser.currentElementNameMatches(
                                                           "MinuitStrategy" ) )
            {
              minuitStrategy = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
            }
            else if( elementParser.currentElementNameMatches(
                                                          "MinuitTolerance" ) )
            {
              minuitToleranceFraction = BOL::StringParser::stringToDouble(
                             elementParser.getTrimmedCurrentElementContent() );
            }
            else if( elementParser.currentElementNameMatches(
                                            "NumberOfPotentialSamplePoints" ) )
            {
              numberOfPotentialSamplePoints = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
            }
            else if( elementParser.currentElementNameMatches(
                                                       "BouncePotentialFit" ) )
            {
              BOL::AsciiXmlParser nestedParser;
              nestedParser.loadString(
                             elementParser.getTrimmedCurrentElementContent() );
              while( nestedParser.readNextElement() )
              {
                if( nestedParser.currentElementNameMatches( "ClassType" ) )
                {
                  pathBouncePotentialFitClass.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
                }
                else if( nestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
                {
                  pathBouncePotentialFitArguments.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
                }
              }
            }
            else if( elementParser.currentElementNameMatches(
                                                     "PathParameterization" ) )
            {
              BOL::AsciiXmlParser nestedParser;
              nestedParser.loadString(
                             elementParser.getTrimmedCurrentElementContent() );
              while( nestedParser.readNextElement() )
              {
                if( nestedParser.currentElementNameMatches( "ClassType" ) )
                {
                  pathParameterizationClass.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
                }
                else if( nestedParser.currentElementNameMatches(
                                                     "ConstructorArguments" ) )
                {
                  pathParameterizationArguments.assign(
                              nestedParser.getTrimmedCurrentElementContent() );
                }
              }
            }
            else if( elementParser.currentElementNameMatches(
                                                    "NodeMovementThreshold" ) )
            {
              nodeMovementThreshold = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
            }
          }

          if( ( pathParameterizationClass.compare( "PathFromNodes" ) == 0 )
              ||
              ( pathParameterizationClass.compare(
                                   "PathFromPolynomialCoefficients" ) == 0 ) )
          {
            size_t numberOfVaryingCoefficientsPerField( 3 );
            size_t numberOfVaryingNodes( 7 );
            std::string nodeParameterizationStyle( "NodesOnParallelPlanes" );
            std::string interpolationStyle( "LinearSplinePath" );
            elementParser.loadString( pathParameterizationArguments );
            while( elementParser.readNextElement() )
            {
              if( elementParser.currentElementNameMatches(
                          "NumberOfCoefficientsPerVaryingFieldBeyondLinear" ) )
              {
                numberOfVaryingCoefficientsPerField
                = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
              }
              else if( elementParser.currentElementNameMatches(
                                                     "NumberOfVaryingNodes" ) )
              {
                numberOfVaryingNodes = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
              }
              else if( elementParser.currentElementNameMatches(
                                                "NodeParameterizationStyle" ) )
              {
                nodeParameterizationStyle.assign(
                             elementParser.getTrimmedCurrentElementContent() );
              }
              else if( elementParser.currentElementNameMatches(
                                                       "InterpolationStyle" ) )
              {
                interpolationStyle.assign(
                             elementParser.getTrimmedCurrentElementContent() );
              }
            }

            TunnelPathFactory* tunnelPathFactory( NULL );
            PathFromNodesFactory* pathFromNodesFactory( NULL );
            if( pathParameterizationClass.compare( "PathFromNodes" ) == 0 )
            {
              NodesFromParameterization* nodesFromParameterization( NULL );
              if( nodeParameterizationStyle.compare( "NodesOnParallelPlanes" )
                  == 0 )
              {
                nodesFromParameterization
                = new NodesOnParallelPlanes( numberOfFields,
                                             numberOfVaryingNodes );
              }
              else if( nodeParameterizationStyle.compare(
                                              "NodesOnBisectingPlanes" ) == 0 )
              {
                nodesFromParameterization
                = new NodesOnBisectingPlanes( numberOfFields,
                                              numberOfVaryingNodes );
              }
              else
              {
                std::stringstream errorStream;
                errorStream
                << "<NodeParameterizationStyle> was not a recognized form! The"
                << " only types currently valid are \"NodesOnParallelPlanes\""
                << " and \"NodesOnBisectingPlanes\".";
                throw std::runtime_error( errorStream.str() );
              }
              if( interpolationStyle.compare( "PolynomialPath" ) == 0 )
              {
                pathFromNodesFactory
                = new PolynomialThroughNodesFactory( numberOfFields,
                                                   nodesFromParameterization );
              }
              else if( interpolationStyle.compare( "LinearSplinePath" ) == 0 )
              {
                pathFromNodesFactory
                = new LinearSplineThroughNodesFactory( numberOfFields,
                                                   nodesFromParameterization );
              }
              else if( interpolationStyle.compare( "QuadraticSplinePath" )
                       == 0 )
              {
                pathFromNodesFactory
                = new QuadraticSplineThroughNodesFactory( numberOfFields,
                                                   nodesFromParameterization );
              }
              tunnelPathFactory = pathFromNodesFactory;
            }
            else if( pathParameterizationClass.compare(
                                      "PathFromPolynomialCoefficients" ) == 0 )
            {
              tunnelPathFactory
              = new PolynomialsFromCoefficientsFactory( numberOfFields,
                                         numberOfVaryingCoefficientsPerField );
            }
            // Now we have the components to assemble the path finder:
            if( tunnelPathFinderClass.compare(
                                        "MinuitNodePotentialMinimizer" ) == 0 )
            {
              bouncePathFinder
              = new MinuitNodePotentialMinimizer( pathFromNodesFactory,
                                                  *potentialFunction,
                                                  movesPerImprovement,
                                                  minuitStrategy,
                                                  minuitToleranceFraction,
                                                  nodeMovementThreshold );
            }
            else if( tunnelPathFinderClass.compare(
                                        "MinuitPathPotentialMinimizer" ) == 0 )
            {
              bouncePathFinder
              = new MinuitPathPotentialMinimizer( tunnelPathFactory,
                                                  *potentialFunction,
                                                 numberOfPotentialSamplePoints,
                                                  movesPerImprovement,
                                                  minuitStrategy,
                                                  minuitToleranceFraction );
            }
            else if( tunnelPathFinderClass.compare(
                                           "MinuitPathBounceMinimizer" ) == 0 )
            {
              // We need to check for a new way of calculating the bounce,
              // which will be used for improving the path but not to
              // calculate the final bounce action along the final path. The
              // input for the final bounce action calculation is used as a
              // default for any input not given explicitly.
              if( pathBouncePotentialFitClass.compare(
                                              "BubbleShootingOnSpline" ) == 0 )
              {
                elementParser.loadString( pathBouncePotentialFitArguments );
                while( elementParser.readNextElement() )
                {
                  if( elementParser.currentElementNameMatches(
                                             "NumberOfNodesForPotentialFit" ) )
                  {
                    numberOfPotentialSegmentsForBounce
                    = ( 1 + BOL::StringParser::stringToInt(
                           elementParser.getTrimmedCurrentElementContent() ) );
                  }
                  else if( elementParser.currentElementNameMatches(
                                               "NumberShootAttemptsAllowed" ) )
                  {
                    shootAttemptsForBounce = BOL::StringParser::stringToInt(
                             elementParser.getTrimmedCurrentElementContent() );
                  }
                  else if( elementParser.currentElementNameMatches(
                                                         "RadialResolution" ) )
                  {
                    lengthScaleResolutionForBounce
                    = BOL::StringParser::stringToDouble(
                             elementParser.getTrimmedCurrentElementContent() );
                  }
                }
                bouncePathFinder
                = new MinuitPathBounceMinimizer( tunnelPathFactory,
                                new BubbleShootingOnSpline( *potentialFunction,
                                            numberOfPotentialSegmentsForBounce,
                                                lengthScaleResolutionForBounce,
                                                      shootAttemptsForBounce ),
                                                 movesPerImprovement,
                                                 minuitStrategy,
                                                 minuitToleranceFraction );
              }
              else
              {
                std::stringstream errorStream;
                errorStream
                << "<BouncePotentialFit> within"
                << " <TunnelPathFinder><ConstructorArguments> was not a"
                << " recognized form! The only type currently valid is"
                << " \"BubbleShootingOnSpline\".";
                throw std::runtime_error( errorStream.str() );
              }
            }
          }
          else
          {
            std::stringstream errorStream;
            errorStream
            << "<TunnelPathFinder> was not a recognized form! The only types"
            << " currently valid are \"PathFromNodes\" and"
            << " \"PathFromPolynomialCoefficients\".";
            throw std::runtime_error( errorStream.str() );
          }
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
        = new BounceAlongPathWithThreshold( *potentialFunction,
                                            bouncePathFinder,
                                            bounceActionCalculator,
                                            tunnelingStrategy,
                                            survivalProbabilityThreshold,
                                            temperatureAccuracy,
                                            evaporationResolution,
                                            thermalIntegrationResolution );
      }
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
    delete deleterForSlhaManager;
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

    slhaManager->UpdateSlhaData( parameterFilename );
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
