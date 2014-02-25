/*
 * VevaciousPlusPlus.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  VevaciousPlusPlus::VevaciousPlusPlus( BOL::ArgumentParser& argumentParser,
                                        PotentialMinimizer& potentialMinimizer,
                                   TunnelingCalculator& tunnelingCalculator ) :
    runTimer( 10.0 ),
    potentialMinimizer( potentialMinimizer ),
    thresholdDistanceSquaredFraction( pow( BOL::StringParser::stringToDouble(
                           argumentParser.fromTag( "MinimaSeparationThreshold",
                                                   "0.05" ) ),
                                           2 ) ),
    extremumSeparationSquaredThreshold( 100.0 ),
    dsbMinimum(),
    globalMinimum(),
    panicMinimum(),
    dsbIsMetastable( false ),
    numberOfNudges( BOL::StringParser::stringToInt( argumentParser.fromTag(
                                                                "SaddleNudges",
                                                                     "1" ) ) ),
    tunnelingCalculator( tunnelingCalculator ),
    quantumLifetimeThreshold( BOL::StringParser::stringToDouble(
                            argumentParser.fromTag( "QuantumLifetimeThreshold",
                                                    "0.217" ) ) ),
    thermalSurvivalThreshold( BOL::StringParser::stringToDouble(
                            argumentParser.fromTag( "ThermalSurvivalThreshold",
                                                    "0.01" ) ) )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "VevaciousPlusPlus::VevaciousPlusPlus( ... )";
    std::cout << std::endl;/**/
  }

  VevaciousPlusPlus::~VevaciousPlusPlus()
  {
    // This does nothing.
  }


  void VevaciousPlusPlus::RunPoint( std::string const& parameterFilename )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "VevaciousPlusPlus::RunPoint( \"" << parameterFilename << "\" )";
    std::cout << std::endl;/**/
  }

  void VevaciousPlusPlus::WriteXmlResults( std::string const& xmlFilename )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "VevaciousPlusPlus::WriteXmlResults( \"" << xmlFilename << "\" )";
    std::cout << std::endl;/**/
  }

  void VevaciousPlusPlus::WriteSlhaResults( std::string const& slhaFilename,
                                            bool const writeWarnings )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "VevaciousPlusPlus::WriteSlhaResults( \"" << slhaFilename
    << ", " << writeWarnings << "\" )";
    std::cout << std::endl;/**/

  }

} /* namespace VevaciousPlusPlus */
