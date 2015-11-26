/*
 * LagrangianParameterManager.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/LagrangianParameterManager.hpp"

namespace VevaciousPlusPlus
{

  LagrangianParameterManager::LagrangianParameterManager() :
    LHPC::BasicObserved()
  {
    // This constructor is just an initialization list.
  }

  LagrangianParameterManager::~LagrangianParameterManager()
  {
    // This does nothing.
  }

  // This returns "everything up to the '['" paired with "everything within the
  // '[' to the ']' (not including the brackets themselves), throwing
  // exceptions if there anything other than either no square brackets, or a
  // single '[' with a single ']' coming after it.
  std::pair< std::string, std::string >
  LagrangianParameterManager::SeparateIndexBracket(
                                    std::string const& stringToSeparate ) const
  {
    size_t openBracket( stringToSeparate.find( '[' ) );
    if( openBracket == std::string::npos )
    {
      return std::make_pair( stringToSeparate,
                             std::string( "" ) );
    }
    size_t closeBracket( stringToSeparate.find_last_of( '[' ) );
    if( ( closeBracket == std::string::npos )
        ||
        ( stringToSeparate.find( openBracket,
                                 '[' ) != std::string::npos )
        ||
        ( stringToSeparate.find_last_of( closeBracket,
                                         '[' ) != std::string::npos ) )
    {
      throw std::runtime_error(
                           "In parsing variable, [...] not closed properly." );
    }
    return std::make_pair( stringToSeparate.substr( 0,
                                                    openBracket ),
                           stringToSeparate.substr( ( openBracket + 1 ),
                             ( stringToSeparate.size() - openBracket - 2 ) ) );
  }

} /* namespace VevaciousPlusPlus */
