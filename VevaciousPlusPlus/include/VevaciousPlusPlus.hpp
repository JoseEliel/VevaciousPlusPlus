/*
 * VevaciousPlusPlus.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef VEVACIOUSPLUSPLUS_HPP_
#define VEVACIOUSPLUSPLUS_HPP_

#include "CommonIncludes.hpp"
#include "VersionInformation.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimizer.hpp"
#include "PotentialMinimization/GradientFromStartingPoints.hpp"
#include "PotentialMinimization/StartingPointFinder.hpp"
#include "PotentialMinimization/HomotopyContinuation/Hom4ps2Runner.hpp"
#include "PotentialMinimization/GradientMinimizer.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitPotentialMinimizer.hpp"
#include "TunnelingCalculation/TunnelingCalculator.hpp"
#include "TunnelingCalculation/BounceActionTunneling/CosmoTransitionsRunner.hpp"
#include "TunnelingCalculation/BounceActionTunneling/BounceAlongPathWithThreshold.hpp"
#include "BounceActionEvaluation/BouncePathFinder.hpp"
#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnPathNormalInertialPotential.hpp"
#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialOnParallelPlanes.hpp"
#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialPerpendicularToPath.hpp"
#include "BounceActionEvaluation/BounceActionCalculator.hpp"
#include "BounceActionEvaluation/BubbleShootingOnPathInFieldSpace.hpp"
#include "LagrangianParameterManagement/RunningParameterManager.hpp"
#include "LagrangianParameterManagement/LagrangianParameterManager.hpp"
#include "LagrangianParameterManagement/LesHouchesAccordBlockEntryManager.hpp"
#include "LagrangianParameterManagement/SlhaBlocksWithSpecialCasesManager.hpp"
#include "LagrangianParameterManagement/SlhaCompatibleWithSarahManager.hpp"
#include "PotentialEvaluation/PotentialFunctions/OldFixedScaleOneLoopPotential.hpp"
#include "PotentialEvaluation/PotentialFunctions/OldPotentialFromPolynomialAndMasses.hpp"
#include "PotentialEvaluation/PotentialFunctions/OldRgeImprovedOneLoopPotential.hpp"


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
    VevaciousPlusPlus( SlhaManager& slhaManager,
                       PotentialMinimizer& potentialMinimizer,
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


    // This runs the point parameterized in the given file.
    void RunPoint( std::string const& parameterFilename );

    // This writes the results as an XML file.
    void WriteXmlResults( std::string const& xmlFilename );

    // This writes the results as an SLHA file.
    void WriteSlhaResults( std::string const& slhaFilename,
                           bool const writeWarnings = true );


  protected:
    // This puts the content of the current element of xmlParser into
    // contentDestination, trimmed of leading and trailing whitespace, if the
    // element's name matches elementName.
    static void
    InterpretElementIfNameMatches( BOL::AsciiXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   std::string& contentDestination );

    // This puts the content of the current element of xmlParser into
    // contentDestination, interpreted as a double represented in ASCII, if the
    // element's name matches elementName.
    static void
    InterpretElementIfNameMatches( BOL::AsciiXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   double& contentDestination );

    // This puts the content of the current element of xmlParser into
    // contentDestination, interpreted as an int represented in ASCII, if the
    // element's name matches elementName.
    static void
    InterpretElementIfNameMatches( BOL::AsciiXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   int& contentDestination );

    // This puts the content of the current element of xmlParser into
    // contentDestination, interpreted as a bool represented by
    // case-insensitive "yes/no" or "y/n" or "true/false" or "t/f" or "0/1",
    // if the element's name matches elementName. If the element content
    // doesn't match any valid input, contentDestination is left untouched.
    static void
    InterpretElementIfNameMatches( BOL::AsciiXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   bool& contentDestination );


    // This reads the current element of outerParser and if its name matches
    // elementName, it puts the contents of the child element <ClassType> into
    // className and <ConstructorArguments> into constructorArguments, both
    // stripped of whitespace.
    static void ReadClassAndArguments( BOL::AsciiXmlParser const& outerParser,
                                       std::string const& elementName,
                                       std::string& className,
                                       std::string& constructorArguments );

    // This interprets the given string as the appropriate element of the
    // TunnelingCalculator::TunnelingStrategy enum.
    static TunnelingCalculator::TunnelingStrategy
    InterpretTunnelingStrategy( std::string& tunnelingStrategy );

    // This throws an exception if the survival probability threshold was
    // outside the range for a valid probability.
    static void
    CheckSurvivalProbabilityThreshold( double survivalProbabilityThreshold );


    SlhaManager* slhaManager;
    RunningParameterManager* ownedSlhaManager;
    PotentialFunction* potentialFunction;
    OldPotentialFromPolynomialAndMasses* ownedPotentialFunction;
    PotentialMinimizer* potentialMinimizer;
    PotentialMinimizer* ownedPotentialMinimizer;
    TunnelingCalculator* tunnelingCalculator;
    TunnelingCalculator* ownedTunnelingCalculator;
    time_t currentTime;


    // This decides on the derived class to use for ownedPotentialFunction and
    // constructs it with the arguments parsed from constructorArguments.
    PotentialFunction* SetUpPotentialFunction( std::string const& className,
                                     std::string const& constructorArguments );

    // This decides on the derived class to use for ownedPotentialMinimizer and
    // constructs it with the arguments parsed from constructorArguments.
    PotentialMinimizer* SetUpPotentialMinimizer( std::string const& className,
                                     std::string const& constructorArguments );

    // This decides on the derived class (based on GradientMinimizer) to use
    // to minimize the potential and constructs it with the arguments parsed
    // from constructorArguments.
    PotentialMinimizer*
    SetUpGradientFromStartingPoints( std::string const& constructorArguments );

    // This parses arguments from constructorArguments and uses them to
    // construct a Hom4ps2Runner instance to use to find starting points.
    StartingPointFinder*
    SetUpHom4ps2Runner( std::string const& constructorArguments );

    // This parses arguments from constructorArguments and uses them to
    // construct a MinuitPotentialMinimizer instance to use to minimize the
    // potential from a set of starting points.
    GradientMinimizer*
    SetUpMinuitPotentialMinimizer( std::string const& constructorArguments );

    // This decides on the derived class to use for ownedTunnelingCalculator
    // and constructs it with the arguments parsed from constructorArguments.
    TunnelingCalculator*
    SetUpTunnelingCalculator( std::string const& className,
                              std::string const& constructorArguments );

    // This parses arguments from constructorArguments and uses them to
    // construct a CosmoTransitionsRunner instance to use to calculate the
    // tunneling from the DSB vacuum to the panic vacuum, if possible.
    TunnelingCalculator*
    SetUpCosmoTransitionsRunner( std::string const& constructorArguments );

    // This parses arguments from constructorArguments and uses them to
    // construct a BounceAlongPathWithThreshold instance to use to calculate
    // the tunneling from the DSB vacuum to the panic vacuum, if possible.
    TunnelingCalculator* SetUpBounceAlongPathWithThreshold(
                                     std::string const& constructorArguments );

    // This parses the XMl of tunnelPathFinders to construct a set of
    // BouncePathFinder instances, filling pathFinders with pointers to them.
    void SetUpBouncePathFinders( std::string const& tunnelPathFinders,
                               std::vector< BouncePathFinder* >& pathFinders );

    // This parses arguments from constructorArguments and uses them to
    // construct a MinuitOnPotentialOnParallelPlanes instance to use to try to
    // extremize the bounce action.
    BouncePathFinder* CreateMinuitOnPotentialOnParallelPlanes(
                                     std::string const& constructorArguments );

    // This parses arguments from constructorArguments and uses them to
    // construct a MinuitOnPotentialPerpendicularToPath instance to use to try
    // to extremize the bounce action.
    BouncePathFinder* CreateMinuitOnPotentialPerpendicularToPath(
                                     std::string const& constructorArguments );

    // This parses arguments from constructorArguments and uses them to
    // construct a MinuitOnPathNormalInertialPotential instance to use to try
    // to extremize the bounce action.
    BouncePathFinder* CreateMinuitOnPathNormalInertialPotential(
                                     std::string const& constructorArguments );

    // This decides on the derived class (based on BounceActionCalculator) to use
    // to extremize the bounce action to use to calculate the tunneling
    // probability and constructs it with the arguments parsed from
    // constructorArguments.
    BounceActionCalculator*
    SetUpBounceActionCalculator( std::string const& className,
                                 std::string const& constructorArguments );

    // This parses arguments from constructorArguments and uses them to
    // construct a BubbleShootingOnPathInFieldSpace instance to use to calculate
    // the bounce action on given paths.
    BounceActionCalculator* SetUpBubbleShootingOnPathInFieldSpace(
                                     std::string const& constructorArguments );
  };




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

  // This puts the content of the current element of xmlParser into
  // contentDestination, interpreted as an int represented in ASCII, if the
  // element's name matches elementName.
  inline void VevaciousPlusPlus::InterpretElementIfNameMatches(
                                          BOL::AsciiXmlParser const& xmlParser,
                                                std::string const& elementName,
                                                      int& contentDestination )
  {
    if( xmlParser.currentElementNameMatches( elementName ) )
    {
      contentDestination = BOL::StringParser::stringToInt(
                             ( xmlParser.getTrimmedCurrentElementContent() ) );
    }
  }

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

  // This throws an exception if the survival probability threshold was outside
  // the range for a valid probability.
  inline void VevaciousPlusPlus::CheckSurvivalProbabilityThreshold(
                                          double survivalProbabilityThreshold )
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

  // This decides on the derived class to use for ownedPotentialMinimizer and
  // constructs it with the arguments parsed from constructorArguments.
  inline PotentialMinimizer*
  VevaciousPlusPlus::SetUpPotentialMinimizer( std::string const& className,
                                      std::string const& constructorArguments )
  {
    if( className.compare( "GradientFromStartingPoints" ) == 0 )
    {
      ownedPotentialMinimizer
      = SetUpGradientFromStartingPoints( constructorArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<MinimizerClass> was not a recognized form! The only type"
      << " currently valid is \"GradientFromStartingPoints\".";
      throw std::runtime_error( errorStream.str() );
    }
    return ownedPotentialMinimizer;
  }

  // This parses arguments from constructorArguments and uses them to construct
  // a Hom4ps2Runner instance to use to find starting points.
  inline StartingPointFinder* VevaciousPlusPlus::SetUpHom4ps2Runner(
                                      std::string const& constructorArguments )
  {
    // The <ConstructorArguments> for this class should have child elements
    // <PathToHom4ps2> and <Hom4ps2Argument>.
    std::string pathToHom4ps2( "./HOM4PS2/" );
    std::string homotopyType( "2" );
    BOL::AsciiXmlParser xmlParser;
    xmlParser.loadString( constructorArguments );
    while( xmlParser.readNextElement() )
    {
      InterpretElementIfNameMatches( xmlParser,
                                     "PathToHom4ps2",
                                     pathToHom4ps2 );
      InterpretElementIfNameMatches( xmlParser,
                                     "Hom4ps2Argument",
                                     homotopyType );
    }
    return new Hom4ps2Runner(
                 *(ownedPotentialFunction->HomotopyContinuationTargetSystem()),
                              pathToHom4ps2,
                              homotopyType );
  }

  // This decides on the derived class to use for ownedTunnelingCalculator and
  // constructs it with the arguments parsed from constructorArguments.
  inline TunnelingCalculator* VevaciousPlusPlus::SetUpTunnelingCalculator(
                                                  std::string const& className,
                                      std::string const& constructorArguments )
  {
    if( className.compare( "CosmoTransitionsRunner" ) == 0 )
    {
      ownedTunnelingCalculator
      = SetUpCosmoTransitionsRunner( constructorArguments );
    }
    else if( className.compare( "BounceAlongPathWithThreshold" ) == 0 )
    {
      ownedTunnelingCalculator
      = SetUpBounceAlongPathWithThreshold( constructorArguments );
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
    return ownedTunnelingCalculator;
  }

  // This decides on the derived class (based on BounceActionCalculator) to use
  // to extremize the bounce action to use to calculate the tunneling
  // probability and constructs it with the arguments parsed from
  // constructorArguments.
  inline BounceActionCalculator*
  VevaciousPlusPlus::SetUpBounceActionCalculator( std::string const& className,
                                      std::string const& constructorArguments )
  {
    if( className.compare( "BubbleShootingOnPathInFieldSpace" ) == 0 )
    {
      return SetUpBubbleShootingOnPathInFieldSpace( constructorArguments );
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<BouncePotentialFit> was not a recognized form! The only type"
      << " currently valid is \"BubbleShootingOnPathInFieldSpace\".";
      throw std::runtime_error( errorStream.str() );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* VEVACIOUSPLUSPLUS_HPP_ */
