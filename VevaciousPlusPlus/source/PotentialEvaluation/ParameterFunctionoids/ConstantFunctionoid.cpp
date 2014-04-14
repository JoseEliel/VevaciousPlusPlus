/*
 * ConstantFunctionoid.cpp
 *
 *  Created on: Mar 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  ConstantFunctionoid::ConstantFunctionoid( double const constantValue,
                                            std::string const& creationString,
                                     std::string const& pythonParameterName ) :
    ParameterFunctionoid( creationString,
                          pythonParameterName,
                          constantValue )
  {
    // This constructor is just an initialization list.
  }

  ConstantFunctionoid::~ConstantFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
