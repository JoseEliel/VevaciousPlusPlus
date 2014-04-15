/*
 * VevaciousPlusPlus.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{
  std::string const VevaciousPlusPlus::versionString( "2.prealpha.0012" );
  std::string const VevaciousPlusPlus::citationString( "[none as yet]" );

  VevaciousPlusPlus::VevaciousPlusPlus( BOL::ArgumentParser& argumentParser,
                                        SlhaManager& slhaManager,
                                        PotentialMinimizer& potentialMinimizer,
                                   TunnelingCalculator& tunnelingCalculator ) :
    // runTimer( 600.0 ),
    slhaManager( slhaManager ),
    potentialMinimizer( potentialMinimizer ),
    tunnelingCalculator( tunnelingCalculator )
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

    slhaManager.UpdateSlhaData( parameterFilename );
    potentialMinimizer.FindMinima( 0.0 );
    if( potentialMinimizer.DsbVacuumIsMetaStable() )
    {
      tunnelingCalculator.CalculateTunneling( potentialMinimizer.DsbVacuum(),
                                            potentialMinimizer.PanicVacuum() );
    }
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
