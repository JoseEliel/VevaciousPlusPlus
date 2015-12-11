/*
 * VevaciousPlusPlus.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef VEVACIOUSPLUSPLUS_HPP_
#define VEVACIOUSPLUSPLUS_HPP_


#include "LagrangianParameterManagement/SlhaCompatibleWithSarahManager.hpp"
#include "LagrangianParameterManagement/SlhaBlocksWithSpecialCasesManager.hpp"
#include "LagrangianParameterManagement/LesHouchesAccordBlockEntryManager.hpp"


#include "LagrangianParameterManagement/LagrangianParameterManager.hpp"
#include "PotentialEvaluation/PotentialFunctions/PotentialFromPolynomialWithMasses.hpp"
#include "PotentialMinimization/PotentialMinimizer.hpp"
#include "TunnelingCalculation/TunnelingCalculator.hpp"
#include "VersionInformation.hpp"
#include "PotentialEvaluation/PotentialFunctions/FixedScaleOneLoopPotential.hpp"
#include "PotentialEvaluation/PotentialFunctions/RgeImprovedOneLoopPotential.hpp"
#include "PotentialMinimization/GradientFromStartingPoints.hpp"
#include "PotentialMinimization/StartingPointFinder.hpp"
#include "PotentialMinimization/GradientMinimizer.hpp"
#include "PotentialMinimization/StartingPointGeneration/PolynomialAtFixedScalesSolver.hpp"
#include "PotentialMinimization/StartingPointGeneration/PolynomialSystemSolver.hpp"
#include "PotentialMinimization/HomotopyContinuation/Hom4ps2Runner.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitPotentialMinimizer.hpp"
#include "TunnelingCalculation/BounceActionTunneling/CosmoTransitionsRunner.hpp"
#include "TunnelingCalculation/BounceActionTunneling/BounceAlongPathWithThreshold.hpp"
#include "BounceActionEvaluation/BouncePathFinder.hpp"
#include "BounceActionEvaluation/BounceActionCalculator.hpp"
#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialOnParallelPlanes.hpp"
#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialPerpendicularToPath.hpp"
#include "BounceActionEvaluation/BubbleShootingOnPathInFieldSpace.hpp"


namespace VevaciousPlusPlus
{
  class VevaciousPlusPlus
  {
  public:
    // This is the constructor for those who know what they are doing, to allow
    // the main function of the program to decide the components and pass them
    // in to the constructor, allowing for custom components without having to
    // edit the VevaciousPlusPlus files. Those wishing to just use the default
    // possibilities can use the other constructor, which takes the name of an
    // initialization file, and then creates the components based on the data
    // in that file.
    VevaciousPlusPlus( PotentialMinimizer& potentialMinimizer,
                       TunnelingCalculator& tunnelingCalculator );

    // This is the constructor that we expect to be used in normal use: it
    // reads in an initialization file in XML with name given by
    // initializationFileName and then assembles the appropriate components for
    // the objects. The body of the constructor is just a lot of statements
    // reading in XML elements and creating new instances of components. It'd
    // be great to be able to use std::unique_ptrs but we're sticking to
    // allowing non-C++11-compliant compilers.
    VevaciousPlusPlus( std::string const& initializationFileName );

    virtual ~VevaciousPlusPlus();


    // This runs the point parameterized by newInput, which for the default
    // case gives the name of a file with the input parameters, but could in
    // principle itself contain all the necessary parameters.
    void RunPoint( std::string const& newInput );

    // This writes the results as an XML file.
    void WriteXmlResults( std::string const& xmlFilename );

    // This writes the results as an SLHA file.
    void WriteLhaResults( std::string const& slhaFilename,
                           bool const writeWarnings = true );


  protected:
    typedef LesHouchesAccordBlockEntryManager LhaParameterManager;
    typedef PotentialFromPolynomialWithMasses OneLoopPotential;
    typedef std::pair< LhaParameterManager*, OneLoopPotential* >
            FullPotentialDescription;

    // This creates a new LagrangianParameterManager and a new
    // PotentialFunction according to the XML elements in the file given by
    // potentialFunctionInitializationFilename and returns pointers to them.
    static FullPotentialDescription CreateFullPotentialDescription(
                  std::string const& potentialFunctionInitializationFilename );

    // This reads the current element of outerParser and if its name matches
    // elementName, it puts the contents of the child element <ClassType> into
    // className and <ConstructorArguments> into constructorArguments, both
    // stripped of whitespace.
    static void ReadClassAndArguments( BOL::AsciiXmlParser const& outerParser,
                                       std::string const& elementName,
                                       std::string& className,
                                       std::string& constructorArguments );

    // This puts the content of the current element of xmlParser into
    // contentDestination, trimmed of leading and trailing whitespace, if the
    // element's name matches elementName.
    static void
    InterpretElementIfNameMatches( BOL::AsciiXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   std::string& contentDestination );

    // This creates a new LagrangianParameterManager based on the given
    // arguments and returns a pointer to it.
    static LesHouchesAccordBlockEntryManager*
    CreateLagrangianParameterManager( std::string const& classChoice,
                                     std::string const& constructorArguments );

    // This creates a new PotentialFunction based on the given arguments and
    // returns a pointer to it.
    static PotentialFromPolynomialWithMasses*
    CreatePotentialFunction( std::string const& classChoice,
                             std::string const& constructorArguments,
                      LagrangianParameterManager& lagrangianParameterManager );

    // This puts the content of the current element of xmlParser into
    // contentDestination, interpreted as a double represented in ASCII, if the
    // element's name matches elementName.
    static void
    InterpretElementIfNameMatches( BOL::AsciiXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   double& contentDestination );

    // This creates a PotentialMinimizer according to the XML elements in the
    // file given by potentialMinimizerInitializationFilename and returns
    // a pointer to it.
    static PotentialMinimizer* CreatePotentialMinimizer(
                          PotentialFromPolynomialWithMasses& potentialFunction,
                 std::string const& potentialMinimizerInitializationFilename );

    // This creates a new PotentialMinimizer based on the given arguments and
    // returns a pointer to it.
    static PotentialMinimizer* CreatePotentialMinimizer(
                         PotentialFromPolynomialWithMasses& potentialFunction,
                                                std::string const& classChoice,
                                     std::string const& constructorArguments );

    // This creates a new GradientFromStartingPoints based on the given
    // arguments and returns a pointer to it.
    static GradientFromStartingPoints* CreateGradientFromStartingPoints(
                          PotentialFromPolynomialWithMasses& potentialFunction,
                                     std::string const& constructorArguments );

    // This creates a new StartingPointFinder based on the given arguments and
    // returns a pointer to it.
    static StartingPointFinder* CreateStartingPointFinder(
                    PotentialFromPolynomialWithMasses const& potentialFunction,
                                                std::string const& classChoice,
                                     std::string const& constructorArguments );

    // This creates a new PolynomialAtFixedScalesSolver based on the given
    // arguments and returns a pointer to it.
    static PolynomialAtFixedScalesSolver* CreatePolynomialAtFixedScalesSolver(
                    PotentialFromPolynomialWithMasses const& potentialFunction,
                                     std::string const& constructorArguments );

    // This puts the content of the current element of xmlParser into
    // contentDestination, interpreted as an int represented in ASCII, if the
    // element's name matches elementName.
    static void
    InterpretElementIfNameMatches( BOL::AsciiXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   unsigned int& contentDestination );

    // This puts the content of the current element of xmlParser into
    // contentDestination, interpreted as a bool represented by
    // case-insensitive "yes/no" or "y/n" or "true/false" or "t/f" or "0/1",
    // if the element's name matches elementName. If the element content
    // doesn't match any valid input, contentDestination is left untouched.
    static void
    InterpretElementIfNameMatches( BOL::AsciiXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   bool& contentDestination );

    // This creates a new PolynomialSystemSolver based on the given arguments
    // and returns a pointer to it.
    static PolynomialSystemSolver*
    CreatePolynomialSystemSolver( std::string const& classChoice,
                                  std::string const& constructorArguments );

    // This creates a new Hom4ps2Runner based on the given arguments and
    // returns a pointer to it.
    static Hom4ps2Runner*
    CreateHom4ps2Runner( std::string const& constructorArguments );

    // This creates a new GradientMinimizer based on the given arguments and
    // returns a pointer to it.
    static GradientMinimizer*
    CreateGradientMinimizer( PotentialFunction const& potentialFunction,
                             std::string const& classChoice,
                             std::string const& constructorArguments );

    // This creates a new MinuitPotentialMinimizer based on the given arguments
    // and returns a pointer to it.
    static MinuitPotentialMinimizer*
    CreateMinuitPotentialMinimizer( PotentialFunction const& potentialFunction,
                                    std::string const& constructorArguments );

    // This creates a TunnelingCalculator according to the XML elements in the
    // file given by tunnelingCalculatorInitializationFilename and returns
    // a pointer to it.
    static TunnelingCalculator* CreateTunnelingCalculator(
                std::string const& tunnelingCalculatorInitializationFilename );

    // This creates a new PotentialMinimizer based on the given arguments and
    // returns a pointer to it.
    static TunnelingCalculator*
    CreateTunnelingCalculator( std::string const& classChoice,
                               std::string const& constructorArguments );

    // This creates a new CosmoTransitionsRunner based on the given arguments
    // and returns a pointer to it.
    static CosmoTransitionsRunner*
    CreateCosmoTransitionsRunner( std::string const& constructorArguments );

    // This throws an exception if the survival probability threshold was
    // outside the range for a valid probability.
    static void CheckSurvivalProbabilityThreshold(
                                   double const survivalProbabilityThreshold );

    // This interprets the given string as the appropriate element of the
    // TunnelingCalculator::TunnelingStrategy enum.
    static TunnelingCalculator::TunnelingStrategy
    InterpretTunnelingStrategy( std::string& tunnelingStrategy );

    // This creates a new CosmoTransitionsRunner based on the given arguments
    // and returns a pointer to it.
    static BounceAlongPathWithThreshold* CreateBounceAlongPathWithThreshold(
                                     std::string const& constructorArguments );

    // This parses the XMl of tunnelPathFinders to construct a set of
    // BouncePathFinder instances, filling pathFinders with pointers to them.
    static std::vector< BouncePathFinder* >
    CreateBouncePathFinders( std::string const& tunnelPathFinders );

    // This parses arguments from constructorArguments and uses them to
    // construct a MinuitOnPotentialOnParallelPlanes instance to use to try to
    // extremize the bounce action.
    static MinuitOnPotentialOnParallelPlanes*
    CreateMinuitOnPotentialOnParallelPlanes(
                                     std::string const& constructorArguments );

    // This parses arguments from constructorArguments and uses them to
    // construct a MinuitOnPotentialPerpendicularToPath instance to use to try
    // to extremize the bounce action.
    static MinuitOnPotentialPerpendicularToPath*
    CreateMinuitOnPotentialPerpendicularToPath(
                                     std::string const& constructorArguments );

    // This creates a new BounceActionCalculator based on the given arguments
    // and returns a pointer to it.
    static BounceActionCalculator*
    CreateBounceActionCalculator( std::string const& classChoice,
                                  std::string const& constructorArguments );

    // This creates a new BubbleShootingOnPathInFieldSpace based on the given
    // arguments and returns a pointer to it.
    static BubbleShootingOnPathInFieldSpace*
    CreateBubbleShootingOnPathInFieldSpace(
                                     std::string const& constructorArguments );


    LagrangianParameterManager* lagrangianParameterManager;
    LesHouchesAccordBlockEntryManager* ownedLagrangianParameterManager;
    PotentialFromPolynomialWithMasses* ownedPotentialFunction;
    PotentialMinimizer* potentialMinimizer;
    PotentialMinimizer* ownedPotentialMinimizer;
    TunnelingCalculator* tunnelingCalculator;
    TunnelingCalculator* ownedTunnelingCalculator;
    time_t currentTime;
  };





  // This reads the current element of outerParser and if its name matches
  // elementName, it puts the contents of the child element <ClassType> into
  // className and <ConstructorArguments> into constructorArguments, both
  // stripped of whitespace.
  inline void VevaciousPlusPlus::ReadClassAndArguments(
                                        BOL::AsciiXmlParser const& outerParser,
                                                std::string const& elementName,
                                                        std::string& className,
                                            std::string& constructorArguments )
  {
    if( outerParser.currentElementNameMatches( elementName ) )
    {
      BOL::AsciiXmlParser innerParser;
      innerParser.loadString( outerParser.getTrimmedCurrentElementContent() );
      while( innerParser.readNextElement() )
      {
        InterpretElementIfNameMatches( innerParser,
                                       "ClassType",
                                       className );
        InterpretElementIfNameMatches( innerParser,
                                       "ConstructorArguments",
                                       constructorArguments );
      }
    }
  }

  // This puts the content of the current element of xmlParser into
  // contentDestination, trimmed of leading and trailing whitespace, if the
  // element's name matches elementName.
  inline void VevaciousPlusPlus::InterpretElementIfNameMatches(
                                          BOL::AsciiXmlParser const& xmlParser,
                                                std::string const& elementName,
                                              std::string& contentDestination )
  {
    if( xmlParser.currentElementNameMatches( elementName ) )
    {
      contentDestination.assign( xmlParser.getTrimmedCurrentElementContent() );
    }
  }

  // This puts the content of the current element of xmlParser into
  // contentDestination, interpreted as a double represented in ASCII, if the
  // element's name matches elementName.
  inline void VevaciousPlusPlus::InterpretElementIfNameMatches(
                                          BOL::AsciiXmlParser const& xmlParser,
                                                std::string const& elementName,
                                                   double& contentDestination )
  {
    if( xmlParser.currentElementNameMatches( elementName ) )
    {
      contentDestination = BOL::StringParser::stringToDouble(
                             ( xmlParser.getTrimmedCurrentElementContent() ) );
    }
  }

  // This creates a PotentialMinimizer according to the XML elements in the
  // file given by potentialMinimizerInitializationFilename and returns
  // a pointer to it.
  inline PotentialMinimizer* VevaciousPlusPlus::CreatePotentialMinimizer(
                          PotentialFromPolynomialWithMasses& potentialFunction,
                  std::string const& potentialMinimizerInitializationFilename )
  {
    BOL::AsciiXmlParser xmlParser;
    xmlParser.openRootElementOfFile(
                                    potentialMinimizerInitializationFilename );
    std::string classChoice( "error" );
    std::string constructorArguments( "error" );
    while( xmlParser.readNextElement() )
    {
      ReadClassAndArguments( xmlParser,
                             "PotentialMinimizerClass",
                             classChoice,
                             constructorArguments );
    }
    return CreatePotentialMinimizer( potentialFunction,
                                     classChoice,
                                     constructorArguments );
  }

  // This creates a new PotentialMinimizer based on the given arguments and
  // returns a pointer to it.
  inline PotentialMinimizer* VevaciousPlusPlus::CreatePotentialMinimizer(
                          PotentialFromPolynomialWithMasses& potentialFunction,
                                                std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "GradientFromStartingPoints" )
    {
      return CreateGradientFromStartingPoints( potentialFunction,
                                               constructorArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<PotentialMinimizerClass> was not a recognized class! The only"
      << " option currently valid is \"GradientFromStartingPoints\".";
      throw std::runtime_error( errorStream.str() );
    }
  }

  // This creates a new StartingPointFinder based on the given arguments and
  // returns a pointer to it.
  inline StartingPointFinder* VevaciousPlusPlus::CreateStartingPointFinder(
                    PotentialFromPolynomialWithMasses const& potentialFunction,
                                                std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "PolynomialAtFixedScalesSolver" )
    {
      return CreatePolynomialAtFixedScalesSolver( potentialFunction,
                                                  constructorArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<StartingPointFinderClass> was not a recognized class! The only"
      << " option currently valid is \"PolynomialAtFixedScalesSolver\".";
      throw std::runtime_error( errorStream.str() );
    }
  }

  // This puts the content of the current element of xmlParser into
  // contentDestination, interpreted as an int represented in ASCII, if the
  // element's name matches elementName.
  inline void VevaciousPlusPlus::InterpretElementIfNameMatches(
                                          BOL::AsciiXmlParser const& xmlParser,
                                                std::string const& elementName,
                                             unsigned int& contentDestination )
  {
    if( xmlParser.currentElementNameMatches( elementName ) )
    {
      contentDestination
      = static_cast< unsigned int >( BOL::StringParser::stringToInt(
                           ( xmlParser.getTrimmedCurrentElementContent() ) ) );
    }
  }

  // This creates a new PolynomialSystemSolver based on the given arguments
  // and returns a pointer to it.
  inline PolynomialSystemSolver*
  VevaciousPlusPlus::CreatePolynomialSystemSolver(
                                                std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "Hom4ps2Runner" )
    {
      return CreateHom4ps2Runner( constructorArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<PolynomialSystemSolver> was not a recognized class! The only"
      << " option currently valid is \"Hom4ps2Runner\".";
      throw std::runtime_error( errorStream.str() );
    }
  }

  // This creates a new Hom4ps2Runner based on the given arguments and
  // returns a pointer to it.
  inline Hom4ps2Runner* VevaciousPlusPlus::CreateHom4ps2Runner(
                                      std::string const& constructorArguments )
  {
    BOL::AsciiXmlParser xmlParser;
    xmlParser.loadString( constructorArguments );
    std::string pathToHom4ps2( "error" );
    std::string homotopyType( "error" );
    double resolutionSize( 1.0 );
    while( xmlParser.readNextElement() )
    {
      InterpretElementIfNameMatches( xmlParser,
                                     "PathToHom4ps2",
                                     pathToHom4ps2 );
      InterpretElementIfNameMatches( xmlParser,
                                     "Hom4ps2Argument",
                                     homotopyType );
      InterpretElementIfNameMatches( xmlParser,
                                     "ResolutionSize",
                                     resolutionSize );
    }
    return new Hom4ps2Runner( pathToHom4ps2,
                              homotopyType,
                              resolutionSize );
  }

  // This creates a new GradientMinimizer based on the given arguments and
  // returns a pointer to it.
  inline GradientMinimizer* VevaciousPlusPlus::CreateGradientMinimizer(
                                    PotentialFunction const& potentialFunction,
                                                std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "MinuitPotentialMinimizer" )
    {
      return CreateMinuitPotentialMinimizer( potentialFunction,
                                             constructorArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<GradientMinimizerClass> was not a recognized class! The only"
      << " option currently valid is \"MinuitPotentialMinimizer\".";
      throw std::runtime_error( errorStream.str() );
    }
  }

  // This creates a TunnelingCalculator according to the XML elements in the
  // file given by tunnelingCalculatorInitializationFilename and returns
  // a pointer to it.
  inline TunnelingCalculator* VevaciousPlusPlus::CreateTunnelingCalculator(
                 std::string const& tunnelingCalculatorInitializationFilename )
  {
    BOL::AsciiXmlParser xmlParser;
    xmlParser.openRootElementOfFile(
                                   tunnelingCalculatorInitializationFilename );
    std::string classChoice( "error" );
    std::string constructorArguments( "error" );
    while( xmlParser.readNextElement() )
    {
      ReadClassAndArguments( xmlParser,
                             "TunnelingClass",
                             classChoice,
                             constructorArguments );
    }
    return CreateTunnelingCalculator( classChoice,
                                      constructorArguments );
  }

  // This creates a new PotentialMinimizer based on the given arguments and
  // returns a pointer to it.
  inline TunnelingCalculator*
  VevaciousPlusPlus::CreateTunnelingCalculator( std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "CosmoTransitionsRunner" )
    {
      return CreateCosmoTransitionsRunner( constructorArguments );
    }
    else if( classChoice == "BounceAlongPathWithThreshold" )
    {
      return CreateBounceAlongPathWithThreshold( constructorArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<TunnelingClass> was not a recognized class! The only"
      << " options currently valid are \"BounceAlongPathWithThreshold\" or"
      << " \"CosmoTransitionsRunner\".";
      throw std::runtime_error( errorStream.str() );
    }
  }

  // This throws an exception if the survival probability threshold was
  // outside the range for a valid probability.
  inline void VevaciousPlusPlus::CheckSurvivalProbabilityThreshold(
                                    double const survivalProbabilityThreshold )
  {
    if( !( ( survivalProbabilityThreshold > 0.0 )
           &&
           ( survivalProbabilityThreshold < 1.0 ) ) )
    {
      std::stringstream errorStream;
      errorStream
      << "<SurvivalProbabilityThreshold> must be greater than 0.0 and"
      << " less than 1.0 to be valid! Instead, the XML element gave \""
      << survivalProbabilityThreshold << "\".";
      throw std::runtime_error( errorStream.str() );
    }
  }

  // This interprets the given string as the appropriate element of the
  // TunnelingCalculator::TunnelingStrategy enum.
  inline TunnelingCalculator::TunnelingStrategy
  VevaciousPlusPlus::InterpretTunnelingStrategy(
                                               std::string& tunnelingStrategy )
  {
    BOL::StringParser::transformToLowercase( tunnelingStrategy );
    if( ( tunnelingStrategy == "defaulttunneling" )
        ||
        ( tunnelingStrategy == "thermalthenquantum" ) )
    {
      return TunnelingCalculator::ThermalThenQuantum;
    }
    else if( tunnelingStrategy == "quantumthenthermal" )
    {
      return TunnelingCalculator::QuantumThenThermal;
    }
    else if( tunnelingStrategy == "justthermal" )
    {
      return TunnelingCalculator::JustThermal;
    }
    else if( tunnelingStrategy == "justquantum" )
    {
      return TunnelingCalculator::JustQuantum;
    }
    return TunnelingCalculator::NoTunneling;
  }

  // This parses arguments from constructorArguments and uses them to
  // construct a MinuitOnPotentialOnParallelPlanes instance to use to try to
  // extremize the bounce action.
  inline MinuitOnPotentialOnParallelPlanes*
  VevaciousPlusPlus::CreateMinuitOnPotentialOnParallelPlanes(
                                      std::string const& constructorArguments )
  {
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

  // This creates a new BounceActionCalculator based on the given arguments
  // and returns a pointer to it.
  inline BounceActionCalculator*
  VevaciousPlusPlus::CreateBounceActionCalculator(
                                                std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "BubbleShootingOnPathInFieldSpace" )
    {
      return CreateBubbleShootingOnPathInFieldSpace( constructorArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<BouncePotentialFit> was not a recognized class! The only"
      << " option currently valid is \"BubbleShootingOnPathInFieldSpace\".";
      throw std::runtime_error( errorStream.str() );
    }
  }

  // This creates a new BubbleShootingOnPathInFieldSpace based on the given
  // arguments and returns a pointer to it.
  inline BubbleShootingOnPathInFieldSpace*
  VevaciousPlusPlus::CreateBubbleShootingOnPathInFieldSpace(
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
#endif /* VEVACIOUSPLUSPLUS_HPP_ */
