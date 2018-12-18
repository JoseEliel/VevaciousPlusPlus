#ifndef __boss__VevaciousPlusPlus_vevacious_1_0_hpp__
#define __boss__VevaciousPlusPlus_vevacious_1_0_hpp__

/*
 * VevaciousPlusPlus.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef VEVACIOUSPLUSPLUS_HPP_
#define VEVACIOUSPLUSPLUS_HPP_

#include "PotentialMinimization/PotentialMinimizer.hpp"
#include "TunnelingCalculation/TunnelingCalculator.hpp"
#include <string>
#include <memory>
#include "LagrangianParameterManagement/LesHouchesAccordBlockEntryManager.hpp"
#include "LagrangianParameterManagement/LagrangianParameterManager.hpp"
#include "PotentialEvaluation/PotentialFunctions/PotentialFromPolynomialWithMasses.hpp"
#include <utility>
#include "LHPC/Utilities/RestrictedXmlParser.hpp"
#include "PotentialMinimization/GradientFromStartingPoints.hpp"
#include "PotentialMinimization/StartingPointFinder.hpp"
#include "PotentialMinimization/StartingPointGeneration/PolynomialAtFixedScalesSolver.hpp"
#include "PotentialMinimization/StartingPointGeneration/PolynomialSystemSolver.hpp"
#include "PotentialMinimization/HomotopyContinuation/Hom4ps2Runner.hpp"
#include "PotentialMinimization/HomotopyContinuation/PHCRunner.hpp"
#include "PotentialMinimization/GradientMinimizer.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitPotentialMinimizer.hpp"
#include "TunnelingCalculation/BounceActionTunneling/CosmoTransitionsRunner.hpp"
#include "TunnelingCalculation/BounceActionTunneling/BounceAlongPathWithThreshold.hpp"
#include "BounceActionEvaluation/BouncePathFinder.hpp"
#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialOnParallelPlanes.hpp"
#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialPerpendicularToPath.hpp"
#include "BounceActionEvaluation/BounceActionCalculator.hpp"
#include "BounceActionEvaluation/BubbleShootingOnPathInFieldSpace.hpp"
#include <ctime>
#include <fstream>
#include "VersionInformation.hpp"
#include "LHPC/Utilities/ParsingUtilities.hpp"
#include <sstream>
#include <stdexcept>
#include "Utilities/WarningLogger.hpp"
#include <iostream>
#include <vector>
#include <cstddef>
#include "LagrangianParameterManagement/SlhaCompatibleWithSarahManager.hpp"
#include "LagrangianParameterManagement/SlhaBlocksWithSpecialCasesManager.hpp"
#include "LagrangianParameterManagement/SARAHManager.hpp"
#include "PotentialEvaluation/PotentialFunctions/FixedScaleOneLoopPotential.hpp"
#include "PotentialEvaluation/PotentialFunctions/RgeImprovedOneLoopPotential.hpp"
#define ENUMS_DECLARED
#include "backend_types/vevacious_1_0/abstract_VevaciousPlusPlus.hpp"
#include "gambit/Backends/abstracttypedefs.hpp"
#include "gambit/Backends/wrappertypedefs.hpp"

namespace VevaciousPlusPlus
{
    namespace Utils {
        //fixes the missing make_unique in std=-c++11
        template<typename T, typename... Args>
        std::unique_ptr <T> make_unique(Args &&... args) {
          return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
        }
           }

  class VevaciousPlusPlus : public virtual Abstract_VevaciousPlusPlus
  {
  public:

      VevaciousPlusPlus(const VevaciousPlusPlus&){
        std::cout
                << std::endl
                << " Copy constructor has ben run!!!! ";
        std::cout << std::endl;

      }

      VevaciousPlusPlus& operator=(const VevaciousPlusPlus&){
        std::cout
                << std::endl
                << " Assignment operator has ben run!!!! ";
        std::cout << std::endl;

      }
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

    virtual ~VevaciousPlusPlus(){
      std::cout
              << std::endl
              << " Vevacious object has died! ";
      std::cout << std::endl;

    };


    // This runs the point parameterized by newInput, which for the default
    // case gives the name of a file with the input parameters, but could in
    // principle itself contain all the necessary parameters.
    void RunPoint( std::string const& newInput );
    
    //This reads in a Slha block and passes it over to LagrangianParameterManager updating 
    //the given parameters in the blockset object. The parameter values are given in a vector 
    // and the dimension is given, with 1 is given the block is read as a list, if n is given
    // the vector is read as a n x n matrix in sequential order row by row. 
    
    void ReadLhaBlock( std::string const& uppercaseBlockName,
    				   double const scale, 
    				   std::vector<std::pair<int,double>> const& parameters, 
    				   int const dimension );

    // This writes the results as an XML file.
    void WriteResultsAsXmlFile( std::string const& xmlFilename );
    
    // This gives the results as a string.
    std::string GetResultsAsString();

     // This gives the lifetime in seconds
     double GetLifetimeInSeconds();

     // This gives the upper bound on thermal survival probability
     double GetThermalProbability();

    // This writes the results as an SLHA file.
    void AppendResultsToLhaFile( std::string const& lhaFilename,
                                 bool const writeWarnings = true );


  protected:
    typedef PotentialFromPolynomialWithMasses OneLoopPotential;
    typedef std::pair<std::unique_ptr <LesHouchesAccordBlockEntryManager>, std::unique_ptr <PotentialFromPolynomialWithMasses> >
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
    static void
    ReadClassAndArguments( LHPC::RestrictedXmlParser const& outerParser,
                           std::string const& elementName,
                           std::string& className,
                           std::string& constructorArguments );

    // This puts the content of the current element of xmlParser into
    // contentDestination, trimmed of leading and trailing whitespace, if the
    // element's name matches elementName.
    static void
    InterpretElementIfNameMatches( LHPC::RestrictedXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   std::string& contentDestination );

    // This creates a new LagrangianParameterManager based on the given
    // arguments and returns a pointer to it.
    static std::unique_ptr<LesHouchesAccordBlockEntryManager>
    CreateLagrangianParameterManager( std::string const& classChoice,
                                     std::string const& constructorArguments );

    // This creates a new PotentialFunction based on the given arguments and
    // returns a pointer to it.
    static std::unique_ptr<PotentialFromPolynomialWithMasses>
    CreatePotentialFunction( std::string const& classChoice,
                             std::string const& constructorArguments,
                      LagrangianParameterManager& lagrangianParameterManager );

    // This puts the content of the current element of xmlParser into
    // contentDestination, interpreted as a double represented in ASCII, if the
    // element's name matches elementName.
    static void
    InterpretElementIfNameMatches( LHPC::RestrictedXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   double& contentDestination );

    // This creates a PotentialMinimizer according to the XML elements in the
    // file given by potentialMinimizerInitializationFilename and returns
    // a pointer to it.
    static std::unique_ptr<PotentialMinimizer> CreatePotentialMinimizer(
                          PotentialFromPolynomialWithMasses& potentialFunction,
                 std::string const& potentialMinimizerInitializationFilename );

    // This creates a new PotentialMinimizer based on the given arguments and
    // returns a pointer to it.
    static std::unique_ptr<PotentialMinimizer> CreatePotentialMinimizer(
                         PotentialFromPolynomialWithMasses& potentialFunction,
                                                std::string const& classChoice,
                                     std::string const& constructorArguments );

    // This creates a new GradientFromStartingPoints based on the given
    // arguments and returns a pointer to it.
    static std::unique_ptr<GradientFromStartingPoints> CreateGradientFromStartingPoints(
                          PotentialFromPolynomialWithMasses& potentialFunction,
                                     std::string const& constructorArguments );

    // This creates a new StartingPointFinder based on the given arguments and
    // returns a pointer to it.
    static std::unique_ptr<StartingPointFinder> CreateStartingPointFinder(
                    PotentialFromPolynomialWithMasses const& potentialFunction,
                                                std::string const& classChoice,
                                     std::string const& constructorArguments );

    // This creates a new PolynomialAtFixedScalesSolver based on the given
    // arguments and returns a pointer to it.
    static std::unique_ptr<PolynomialAtFixedScalesSolver> CreatePolynomialAtFixedScalesSolver(
                    PotentialFromPolynomialWithMasses const& potentialFunction,
                                     std::string const& constructorArguments );

    // This puts the content of the current element of xmlParser into
    // contentDestination, interpreted as an int represented in ASCII, if the
    // element's name matches elementName.
    static void
    InterpretElementIfNameMatches( LHPC::RestrictedXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   unsigned int& contentDestination );

    // This puts the content of the current element of xmlParser into
    // contentDestination, interpreted as a bool represented by
    // case-insensitive "yes/no" or "y/n" or "true/false" or "t/f" or "0/1",
    // if the element's name matches elementName. If the element content
    // doesn't match any valid input, contentDestination is left untouched.
    static void
    InterpretElementIfNameMatches( LHPC::RestrictedXmlParser const& xmlParser,
                                   std::string const& elementName,
                                   bool& contentDestination );

    // This creates a new PolynomialSystemSolver based on the given arguments
    // and returns a pointer to it.
    static std::unique_ptr<PolynomialSystemSolver>
    CreatePolynomialSystemSolver( std::string const& classChoice,
                                  std::string const& constructorArguments );

    // This creates a new Hom4ps2Runner based on the given arguments and
    // returns a pointer to it.
    static std::unique_ptr<Hom4ps2Runner>
    CreateHom4ps2Runner( std::string const& constructorArguments );

    // This creates a new PHCRunner based on the given arguments and
    // returns a pointer to it.
    static std::unique_ptr<PHCRunner>
    CreatePHCRunner( std::string const& constructorArguments );

    // This creates a new GradientMinimizer based on the given arguments and
    // returns a pointer to it.
    static std::unique_ptr<GradientMinimizer>
    CreateGradientMinimizer( PotentialFunction const& potentialFunction,
                             std::string const& classChoice,
                             std::string const& constructorArguments );

    // This creates a new MinuitPotentialMinimizer based on the given arguments
    // and returns a pointer to it.
    static std::unique_ptr<MinuitPotentialMinimizer>
    CreateMinuitPotentialMinimizer( PotentialFunction const& potentialFunction,
                                    std::string const& constructorArguments );

    // This creates a TunnelingCalculator according to the XML elements in the
    // file given by tunnelingCalculatorInitializationFilename and returns
    // a pointer to it.
    static std::unique_ptr<TunnelingCalculator> CreateTunnelingCalculator(
                std::string const& tunnelingCalculatorInitializationFilename );

    // This creates a new PotentialMinimizer based on the given arguments and
    // returns a pointer to it.
    static std::unique_ptr<TunnelingCalculator>
    CreateTunnelingCalculator( std::string const& classChoice,
                               std::string const& constructorArguments );

    // This creates a new CosmoTransitionsRunner based on the given arguments
    // and returns a pointer to it.
    static std::unique_ptr<CosmoTransitionsRunner>
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
    static std::unique_ptr<BounceAlongPathWithThreshold> CreateBounceAlongPathWithThreshold(
                                     std::string const& constructorArguments );

    // This parses the XMl of tunnelPathFinders to construct a set of
    // BouncePathFinder instances, filling pathFinders with pointers to them.
    static std::vector< std::unique_ptr<BouncePathFinder> >
    CreateBouncePathFinders( std::string const& tunnelPathFinders );

    // This parses arguments from constructorArguments and uses them to
    // construct a MinuitOnPotentialOnParallelPlanes instance to use to try to
    // extremize the bounce action.
    static std::unique_ptr<MinuitOnPotentialOnParallelPlanes>
     CreateMinuitOnPotentialOnParallelPlanes(
                                     std::string const& constructorArguments );

    // This parses arguments from constructorArguments and uses them to
    // construct a MinuitOnPotentialPerpendicularToPath instance to use to try
    // to extremize the bounce action.
    static std::unique_ptr<MinuitOnPotentialPerpendicularToPath>
    CreateMinuitOnPotentialPerpendicularToPath(
                                     std::string const& constructorArguments );

    // This creates a new BounceActionCalculator based on the given arguments
    // and returns a pointer to it.
    static std::unique_ptr<BounceActionCalculator>
    CreateBounceActionCalculator( std::string const& classChoice,
                                  std::string const& constructorArguments );

    // This creates a new BubbleShootingOnPathInFieldSpace based on the given
    // arguments and returns a pointer to it.
    static std::unique_ptr<BubbleShootingOnPathInFieldSpace>
    CreateBubbleShootingOnPathInFieldSpace(
                                     std::string const& constructorArguments );


    std::unique_ptr<LagrangianParameterManager> lagrangianParameterManager;
    std::unique_ptr<PotentialFromPolynomialWithMasses> ownedPotentialFunction;
    std::unique_ptr<PotentialMinimizer> potentialMinimizer;
    std::unique_ptr<TunnelingCalculator> tunnelingCalculator;
    std::vector< std::string > warningMessagesFromConstructor;
    std::string resultsFromLastRunAsXml;
    std::vector< std::string > warningMessagesFromLastRun;


    // This prepares the results in XML format, stored in resultsAsXml;
    void PrepareResultsAsXml();

    // This returns a vector which is the union of
    // warningMessagesFromConstructor with warningMessagesFromLastRun.
    std::vector< std::string > WarningMessagesToReport() const;

  public:
      Abstract_VevaciousPlusPlus* pointer_copy__BOSS();

      using Abstract_VevaciousPlusPlus::pointer_assign__BOSS;
      void pointer_assign__BOSS(Abstract_VevaciousPlusPlus* in);


  public:
      void AppendResultsToLhaFile__BOSS(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >&);

};





  // This writes the results as an XML file.
  inline void
  VevaciousPlusPlus::WriteResultsAsXmlFile( std::string const& xmlFilename )
  {
    std::time_t currentTime( time( NULL ) );
    std::ofstream xmlFile( xmlFilename.c_str() );
    xmlFile << "<VevaciousResults>\n"
    "  <ReferenceData>\n"
    "     <VevaciousVersion>\n"
    "       " << VersionInformation::CurrentVersion() << "\n"
    "     </VevaciousVersion>\n"
    "     <CitationArticle>\n"
    "       " << VersionInformation::CurrentCitation() << "\n"
    "     </CitationArticle>\n"
    "     <ResultTimestamp>\n"
    "       " << std::string( ctime( &currentTime ) )
    << "     </ResultTimestamp>\n"
    "  </ReferenceData>\n"
    << resultsFromLastRunAsXml << "\n"
    << "</VevaciousResults>\n";
    xmlFile.close();
    std::cout << std::endl << "Wrote results in XML in file \"" << xmlFilename
    << "\"." << std::endl;
  }

  inline std::string VevaciousPlusPlus::GetResultsAsString()
  {
    std::string result= "Error";
    if( potentialMinimizer->DsbVacuumIsStable() )
    {
     result="Stable";
    }
    else
    {
     result="Metastable";
    }
    return result;
  }

     // This gives the Lifetime in seconds as output.
    inline double
    VevaciousPlusPlus::GetLifetimeInSeconds() {
      if (tunnelingCalculator->QuantumSurvivalProbability() >= 0.0)
      {
      return tunnelingCalculator->QuantumLifetimeInSeconds();
      }
    else
      {
        return -1;
      }
    }


    // This gives the upper bound on the thermal survival probability as output.
    inline double
    VevaciousPlusPlus::GetThermalProbability() {
      if (tunnelingCalculator->ThermalSurvivalProbability()  >= 0.0)
      {
        return tunnelingCalculator->ThermalSurvivalProbability();
      }
      else
      {
        return -1;
      }
    }

      // This reads the current element of outerParser and if its name matches
  // elementName, it puts the contents of the child element <ClassType> into
  // className and <ConstructorArguments> into constructorArguments, both
  // stripped of whitespace.
  inline void VevaciousPlusPlus::ReadClassAndArguments(
                                  LHPC::RestrictedXmlParser const& outerParser,
                                                std::string const& elementName,
                                                        std::string& className,
                                            std::string& constructorArguments )
  {
    if( outerParser.CurrentName() == elementName )
    {
      LHPC::RestrictedXmlParser innerParser;
      innerParser.LoadString( outerParser.TrimmedCurrentBody() );
      while( innerParser.ReadNextElement() )
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
                                    LHPC::RestrictedXmlParser const& xmlParser,
                                                std::string const& elementName,
                                              std::string& contentDestination )
  {
    if( xmlParser.CurrentName() == elementName )
    {
      contentDestination.assign( xmlParser.TrimmedCurrentBody() );
    }
  }

  // This puts the content of the current element of xmlParser into
  // contentDestination, interpreted as a double represented in ASCII, if the
  // element's name matches elementName.
  inline void VevaciousPlusPlus::InterpretElementIfNameMatches(
                                    LHPC::RestrictedXmlParser const& xmlParser,
                                                std::string const& elementName,
                                                   double& contentDestination )
  {
    if( xmlParser.CurrentName() == elementName )
    {
      contentDestination = LHPC::ParsingUtilities::StringToDouble(
                                              xmlParser.TrimmedCurrentBody() );
    }
  }

  // This creates a PotentialMinimizer according to the XML elements in the
  // file given by potentialMinimizerInitializationFilename and returns
  // a pointer to it.
  inline std::unique_ptr<PotentialMinimizer> VevaciousPlusPlus::CreatePotentialMinimizer(
                          PotentialFromPolynomialWithMasses& potentialFunction,
                  std::string const& potentialMinimizerInitializationFilename )
  {
    LHPC::RestrictedXmlParser xmlParser;
    xmlParser.OpenRootElementOfFile(
                                    potentialMinimizerInitializationFilename );
    std::string classChoice( "error" );
    std::string constructorArguments( "error" );
    while( xmlParser.ReadNextElement() )
    {
      ReadClassAndArguments( xmlParser,
                             "PotentialMinimizerClass",
                             classChoice,
                             constructorArguments );
    }
    return std::move(CreatePotentialMinimizer( potentialFunction,
                                     classChoice,
                                     constructorArguments ));
  }

  // This creates a new PotentialMinimizer based on the given arguments and
  // returns a pointer to it.
  inline std::unique_ptr<PotentialMinimizer> VevaciousPlusPlus::CreatePotentialMinimizer(
                          PotentialFromPolynomialWithMasses& potentialFunction,
                                                std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "GradientFromStartingPoints" )
    {
      return std::move(CreateGradientFromStartingPoints( potentialFunction,
                                               constructorArguments ));
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
  inline std::unique_ptr<StartingPointFinder> VevaciousPlusPlus::CreateStartingPointFinder(
                    PotentialFromPolynomialWithMasses const& potentialFunction,
                                                std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "PolynomialAtFixedScalesSolver" )
    {
      return std::move(CreatePolynomialAtFixedScalesSolver( potentialFunction,
                                                  constructorArguments ));
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
                                    LHPC::RestrictedXmlParser const& xmlParser,
                                                std::string const& elementName,
                                             unsigned int& contentDestination )
  {
    if( xmlParser.CurrentName() == elementName )
    {
      contentDestination = static_cast< unsigned int >(
                                    LHPC::ParsingUtilities::BaseTenStringToInt(
                                            xmlParser.TrimmedCurrentBody() ) );
    }
  }

  // This creates a new PolynomialSystemSolver based on the given arguments
  // and returns a pointer to it.
  inline std::unique_ptr<PolynomialSystemSolver>
  VevaciousPlusPlus::CreatePolynomialSystemSolver(
                                                std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "Hom4ps2Runner" )
    {
      return std::move(CreateHom4ps2Runner( constructorArguments ));
    }
	else if( classChoice == "PHCRunner" )
    {
      return std::move(CreatePHCRunner( constructorArguments ));
    }
    else
    {
      std::stringstream errorStream;
      errorStream
      << "<PolynomialSystemSolver> was not a recognized class! The only"
      << " options currently valid are \"Hom4ps2Runner\" or \"PHCRunner\"." << std::endl;
	  errorStream << "Classchoice: " << classChoice << std::endl << "Constructorarguments:" << constructorArguments<< std::endl;
      throw std::runtime_error( errorStream.str() );
    }
  }

  // This creates a new Hom4ps2Runner based on the given arguments and
  // returns a pointer to it.
  inline std::unique_ptr<Hom4ps2Runner> VevaciousPlusPlus::CreateHom4ps2Runner(
                                      std::string const& constructorArguments )
  {
    LHPC::RestrictedXmlParser xmlParser;
    xmlParser.LoadString( constructorArguments );
    std::string pathToHom4ps2( "error" );
    std::string homotopyType( "error" );
    double resolutionSize( 1.0 );
    while( xmlParser.ReadNextElement() )
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
    return Utils::make_unique<Hom4ps2Runner>( pathToHom4ps2,
                              homotopyType,
                              resolutionSize );
  }
  
    inline std::unique_ptr<PHCRunner> VevaciousPlusPlus::CreatePHCRunner(
                                      std::string const& constructorArguments )
  {
    LHPC::RestrictedXmlParser xmlParser;
    xmlParser.LoadString( constructorArguments );
	std::string pathToPHC( "error" );
    double resolutionSize( 1.0 );
	unsigned int taskcount(1);
    while( xmlParser.ReadNextElement() )
    {
	InterpretElementIfNameMatches( xmlParser,
                                     "PathToPHC",
                                     pathToPHC );
	InterpretElementIfNameMatches( xmlParser,
                                     "ResolutionSize",
                                     resolutionSize );
	InterpretElementIfNameMatches( xmlParser,
                                     "Tasks",
                                     taskcount );
    }
    return Utils::make_unique<PHCRunner>(  pathToPHC, resolutionSize, taskcount);
  }

  // This creates a new GradientMinimizer based on the given arguments and
  // returns a pointer to it.
  inline std::unique_ptr<GradientMinimizer> VevaciousPlusPlus::CreateGradientMinimizer(
                                    PotentialFunction const& potentialFunction,
                                                std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "MinuitPotentialMinimizer" )
    {
      return std::move(CreateMinuitPotentialMinimizer( potentialFunction,
                                             constructorArguments ));
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
  inline std::unique_ptr<TunnelingCalculator> VevaciousPlusPlus::CreateTunnelingCalculator(
                 std::string const& tunnelingCalculatorInitializationFilename )
  {
    LHPC::RestrictedXmlParser xmlParser;
    xmlParser.OpenRootElementOfFile(
                                   tunnelingCalculatorInitializationFilename );
    std::string classChoice( "error" );
    std::string constructorArguments( "error" );
    while( xmlParser.ReadNextElement() )
    {
      ReadClassAndArguments( xmlParser,
                             "TunnelingClass",
                             classChoice,
                             constructorArguments );
    }
    return std::move(CreateTunnelingCalculator( classChoice,
                                      constructorArguments ));
  }

  // This creates a new PotentialMinimizer based on the given arguments and
  // returns a pointer to it.
  inline std::unique_ptr<TunnelingCalculator>
  VevaciousPlusPlus::CreateTunnelingCalculator( std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "CosmoTransitionsRunner" )
    {
      return std::move(CreateCosmoTransitionsRunner( constructorArguments ));
    }
    else if( classChoice == "BounceAlongPathWithThreshold" )
    {
      return std::move(CreateBounceAlongPathWithThreshold( constructorArguments ));
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

  // This parses arguments from constructorArguments and uses them to
  // construct a MinuitOnPotentialOnParallelPlanes instance to use to try to
  // extremize the bounce action.
  inline std::unique_ptr<MinuitOnPotentialOnParallelPlanes>
  VevaciousPlusPlus::CreateMinuitOnPotentialOnParallelPlanes(
                                      std::string const& constructorArguments )
  {
    unsigned int numberOfPathSegments( 100 );
    unsigned int minuitStrategy( 1 );
    double minuitToleranceFraction( 0.5 );
    LHPC::RestrictedXmlParser xmlParser;
    xmlParser.LoadString( constructorArguments );
    while( xmlParser.ReadNextElement() )
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
    return Utils::make_unique<MinuitOnPotentialOnParallelPlanes>(numberOfPathSegments,
                                                  minuitStrategy,
                                                  minuitToleranceFraction );
  }

  // This creates a new BounceActionCalculator based on the given arguments
  // and returns a pointer to it.
  inline std::unique_ptr<BounceActionCalculator>
  VevaciousPlusPlus::CreateBounceActionCalculator(
                                                std::string const& classChoice,
                                      std::string const& constructorArguments )
  {
    if( classChoice == "BubbleShootingOnPathInFieldSpace" )
    {
      return std::move(CreateBubbleShootingOnPathInFieldSpace( constructorArguments ));
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
  inline std::unique_ptr<BubbleShootingOnPathInFieldSpace>
  VevaciousPlusPlus::CreateBubbleShootingOnPathInFieldSpace(
                                      std::string const& constructorArguments )
  {
    double lengthScaleResolutionForBounce( 0.05 );
    unsigned int shootAttemptsForBounce( 32 );
    LHPC::RestrictedXmlParser xmlParser;
    xmlParser.LoadString( constructorArguments );
    while( xmlParser.ReadNextElement() )
    {
      InterpretElementIfNameMatches( xmlParser,
                                     "RadialResolution",
                                     lengthScaleResolutionForBounce );
      InterpretElementIfNameMatches( xmlParser,
                                     "NumberShootAttemptsAllowed",
                                     shootAttemptsForBounce );
    }
    return
    Utils::make_unique<BubbleShootingOnPathInFieldSpace>( lengthScaleResolutionForBounce,
                                          shootAttemptsForBounce );
  }

  // This returns a vector which is the union of
  // warningMessagesFromConstructor with warningMessagesFromLastRun.
  inline std::vector< std::string >
  VevaciousPlusPlus::WarningMessagesToReport() const
  {
    std::vector< std::string >
    warningMessagesToReport( warningMessagesFromConstructor );
    warningMessagesToReport.insert( warningMessagesToReport.end(),
                                    warningMessagesFromLastRun.begin(),
                                    warningMessagesFromLastRun.end() );
    return warningMessagesToReport;
  }


} /* namespace VevaciousPlusPlus */
#endif /* VEVACIOUSPLUSPLUS_HPP_ */

#endif /* __boss__VevaciousPlusPlus_vevacious_1_0_hpp__ */
