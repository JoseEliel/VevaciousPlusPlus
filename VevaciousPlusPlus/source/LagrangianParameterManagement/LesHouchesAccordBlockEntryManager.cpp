/*
 * LesHouchesAccordBlockEntryManager.cpp
 *
 *  Created on: Oct 29, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/LesHouchesAccordBlockEntryManager.hpp"

namespace VevaciousPlusPlus
{
  std::string const
  LesHouchesAccordBlockEntryManager::blockNameSeparationCharacters(
                                                                 " ,;\t\n\r" );

  LesHouchesAccordBlockEntryManager::LesHouchesAccordBlockEntryManager(
                                       std::string const& validBlocksString ) :
    LagrangianParameterManager(),
    numberOfDistinctActiveParameters( 0 ),
    activeParametersToIndices(),
    activeInterpolatedParameters(),
    validBlocks(),
    lhaParser()
  {
    size_t wordStart( validBlocksString.find_first_of(
                                             blockNameSeparationCharacters ) );
    size_t wordEnd( 0 );
    // If there are any more chars in validBlocks that are not in
    // blockSeparators, we have at least one substring to add.
    while( wordStart != std::string::npos )
    {
      wordEnd = validBlocksString.find_first_of( blockNameSeparationCharacters,
                                                 wordStart );
      validBlocks.insert( validBlocksString.substr( wordStart,
                                                   ( wordEnd - wordStart ) ) );
      if( wordEnd == std::string::npos )
      {
        break;
      }
      wordStart
      = validBlocksString.find_first_not_of( blockNameSeparationCharacters,
                                             wordEnd );
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "LesHouchesAccordBlockEntryManager using SlhaPolynomialFitBlockEntry"
        " objects! Probably should change this to"
        " SlhaLinearlyInterpolatedBlockEntry objects.";
    std::cout << std::endl;/**/
  }

  LesHouchesAccordBlockEntryManager::LesHouchesAccordBlockEntryManager(
                                 std::set< std::string > const& validBlocks ) :
    LagrangianParameterManager(),
    numberOfDistinctActiveParameters( 0 ),
    activeParametersToIndices(),
    activeInterpolatedParameters(),
    validBlocks( validBlocks ),
    lhaParser()
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "LesHouchesAccordBlockEntryManager using SlhaPolynomialFitBlockEntry"
        " objects! Probably should change this to"
        " SlhaLinearlyInterpolatedBlockEntry objects.";
    std::cout << std::endl;/**/
  }

  LesHouchesAccordBlockEntryManager::~LesHouchesAccordBlockEntryManager()
  {
    // This does nothing.
  }


  // This ensures that the given parameter exists in
  // activeInterpolatedParameters and activeParametersToIndices, and returns a
  // reference to its SlhaInterpolatedParameterFunctionoid.
  SlhaInterpolatedParameterFunctionoid const&
  LesHouchesAccordBlockEntryManager::RegisterBlockEntry(
                                             std::string const& parameterName )
  {
    std::map< std::string, size_t >::const_iterator
    alreadyExistsResult( activeParametersToIndices.find( parameterName ) );
    if( alreadyExistsResult != activeParametersToIndices.end() )
    {
      // If the parameter is already active, we unfortunately need to do a
      // linear search through activeInterpolatedParameters looking for the
      // LhaBlockEntryInterpolator with the matching index.
      for( std::vector< LhaBlockEntryInterpolator >::const_iterator
           activeBlockParameter( activeInterpolatedParameters.begin() );
           activeBlockParameter < activeInterpolatedParameters.end();
           ++activeBlockParameter )
      {
        if( activeBlockParameter->IndexInValuesVector()
            == alreadyExistsResult->second )
        {
          return (*activeBlockParameter);
        }
      }
      std::stringstream errorBuilder;
      errorBuilder
      << "LesHouchesAccordBlockEntryManager::RegisterBlockEntry found the"
      << " index " << alreadyExistsResult->second << " for \""
      << parameterName << "\" but no LhaBlockEntryInterpolator in"
      << " activeInterpolatedParameters had that index.";
      throw std::out_of_range( errorBuilder.str() );
    }

    // If the parameter wasn't already active, we add it.
    return AddNewBlockEntry( parameterName );
  }

} /* namespace VevaciousPlusPlus */
