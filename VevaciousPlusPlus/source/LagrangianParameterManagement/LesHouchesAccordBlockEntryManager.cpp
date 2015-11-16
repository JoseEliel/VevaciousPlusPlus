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
                                          std::string const& validBlocksString,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                      std::string const& fixedScaleArgument ) :
    LagrangianParameterManager(),
    numberOfDistinctActiveParameters( 0 ),
    activeParametersToIndices(),
    activeInterpolatedParameters(),
    validBlocks(),
    lhaParser(),
    minimumScaleType( minimumScaleType ),
    minimumScaleArgument( minimumScaleArgument ),
    fixedScaleType( fixedScaleType ),
    fixedScaleArgument( fixedScaleArgument )
  {
    size_t wordStart( validBlocksString.find_first_not_of(
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
                                 std::set< std::string > const& validBlocksSet,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                      std::string const& fixedScaleArgument ) :
    LagrangianParameterManager(),
    numberOfDistinctActiveParameters( 0 ),
    activeParametersToIndices(),
    activeInterpolatedParameters(),
    validBlocks( validBlocksSet ),
    lhaParser(),
    minimumScaleType( minimumScaleType ),
    minimumScaleArgument( minimumScaleArgument ),
    fixedScaleType( fixedScaleType ),
    fixedScaleArgument( fixedScaleArgument )
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


  // This puts all variables with index brackets into a consistent form,
  // putting all those which are a valid block name followed by index brackets
  // into uppercase, to account for SLHA block name case insensitivity.
  std::string LesHouchesAccordBlockEntryManager::FormatVariable(
                                    std::string const& variableToFormat ) const
  {
    std::string const
    trimmedVariable( LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack(
                                                          variableToFormat ) );
    size_t openBracket( trimmedVariable.find( '[' ) );
    if( openBracket == std::string::npos )
    {
      return trimmedVariable;
    }
    std::stringstream formattedStream;
    std::string const blockName( BlockNamePart( trimmedVariable ) );
    if( validBlocks.find( blockName ) != validBlocks.end() )
    {
      formattedStream << blockName;
    }
    else
    {
      formattedStream << trimmedVariable.substr( 0,
                                                 openBracket );
    }
    formattedStream
    << '['
    << FormatIndexBracketContent( trimmedVariable.substr( ( openBracket + 1 ),
                               ( trimmedVariable.size() - openBracket - 2 ) ) )
    << ']';
    return formattedStream.str();
  }

  // This is mainly for debugging.
  std::string LesHouchesAccordBlockEntryManager::AsDebuggingString()
  {
    std::stringstream stringBuilder;

    stringBuilder
    << "numberOfDistinctActiveParameters = "
    << numberOfDistinctActiveParameters << std::endl
    << "activeParametersToIndices = { ";
    for( std::map< std::string, size_t >::const_iterator
         parameterToIndex( activeParametersToIndices.begin() );
         parameterToIndex != activeParametersToIndices.end();
         ++parameterToIndex )
    {
      if( parameterToIndex != activeParametersToIndices.begin() )
      {
        stringBuilder << "," << std::endl;
      }
      stringBuilder
      << "[ \"" << parameterToIndex->first << "\", "
      << parameterToIndex->second << " ]";
    }
    stringBuilder << " }" << std::endl;
    stringBuilder
    << "activeInterpolatedParameters.size() = "
    << activeInterpolatedParameters.size()
    << " (not in the mood to write code to display contents right now)"
    << std::endl
    << "validBlocks = { ";
    for( std::set< std::string >::const_iterator
         validBlock( validBlocks.begin() );
         validBlock != validBlocks.end();
         ++validBlock )
    {
      if( validBlock != validBlocks.begin() )
      {
        stringBuilder << "," << std::endl;
      }
      stringBuilder
      << "\"" << *validBlock << "\"";
    }
    stringBuilder << " }"
    << std::endl << "minimumScaleType = \"" << minimumScaleType << "\""
    << std::endl << "minimumScaleArgument = \"" << minimumScaleArgument << "\""
    << std::endl << "fixedScaleType = \"" << fixedScaleType << "\""
    << std::endl << "fixedScaleArgument = \"" << fixedScaleArgument << "\""
    << std::endl;
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
