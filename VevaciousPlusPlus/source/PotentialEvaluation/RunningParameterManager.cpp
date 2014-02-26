/*
 * RunningParameterManager.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  RunningParameterManager::RunningParameterManager()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RunningParameterManager::RunningParameterManager()";
    std::cout << std::endl;/**/
  }

  RunningParameterManager::~RunningParameterManager()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RunningParameterManager::~RunningParameterManager()";
    std::cout << std::endl;/**/
  }


  // This adds to the list of SLHA blocks that will be read in when updating
  // running parameters. An alias for the block can be provided.
  void
  RunningParameterManager::AddValidSlhaBlock( std::string const& slhaBlock,
                                              std::string const& aliasString )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RunningParameterManager::AddValidSlhaBlock( \"" << slhaBlock
    << "\", \"" << aliasString << "\" )";
    std::cout << std::endl;/**/
  }

  // This creates a new ParameterFunctionoid and a string mapping to it.
  void RunningParameterManager::CreateDerivedParameter(
                                           std::string const& parameterString )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RunningParameterManager::CreateDerivedParameter( \"" << parameterString
    << "\" )";
    std::cout << std::endl;/**/
  }

  // This returns a pointer to the ParameterFunctionoid mapped to by
  // parameterString, returning NULL if there is none.
  ParameterFunctionoid*
  RunningParameterManager::GetFunctionoid( std::string const& parameterString )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RunningParameterManager::GetFunctionoid( \"" << parameterString
    << "\" )";
    std::cout << std::endl;
    return NULL;/**/
  }

} /* namespace VevaciousPlusPlus */
