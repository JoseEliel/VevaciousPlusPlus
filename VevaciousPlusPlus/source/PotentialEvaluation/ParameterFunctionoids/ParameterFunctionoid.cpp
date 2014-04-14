/*
 * ParameterFunctionoid.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  ParameterFunctionoid::ParameterFunctionoid(
                                             std::string const& creationString,
                                        std::string const& pythonParameterName,
                                              double const currentValue ) :
    currentValue( currentValue ),
    creationString( creationString ),
    pythonParameterName( pythonParameterName )
  {
    // This constructor is just an initialization list.
  }

  ParameterFunctionoid::~ParameterFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
