/*
 * SlhaFunctionoid.cpp
 *
 *  Created on: Mar 3, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  SlhaFunctionoid::SlhaFunctionoid( std::string const& indexString ) :
    ParameterFunctionoid(),
    slhaBlock( NULL ),
    indexVector( BOL::StringParser::stringToIntVector( indexString ) ),
    scaleLogarithmPowerCoefficients( 1,
                                     NAN )
  {
    // This constructor is just an initialization list.

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "SlhaFunctionoid::SlhaFunctionoid( \"" << indexString << "\" )"
    << std::endl
    << "slhaBlock = " << slhaBlock << ", indexVector =";
    for( std::vector< int >::const_iterator
         whichIndex( indexVector.begin() );
         whichIndex < indexVector.end();
         ++whichIndex )
    {
      std::cout << " " << *whichIndex;
    }
    std::cout << std::endl;/**/
  }

  SlhaFunctionoid::~SlhaFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
