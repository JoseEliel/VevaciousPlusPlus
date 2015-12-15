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
                                         std::string const& fixedScaleArgument,
                                           std::string const& maximumScaleType,
                                    std::string const& maximumScaleArgument ) :
    LagrangianParameterManager(),
    numberOfDistinctActiveParameters( 0 ),
    activeParametersToIndices(),
    referenceSafeActiveParameters(),
    referenceUnsafeActiveParameters(),
    validBlocks(),
    lhaParser(),
    minimumScaleType( minimumScaleType ),
    minimumScaleArgument( minimumScaleArgument ),
    fixedScaleType( fixedScaleType ),
    fixedScaleArgument( fixedScaleArgument ),
    maximumScaleType( maximumScaleType ),
    maximumScaleArgument( maximumScaleArgument )
  {
    ParseValidBlocks( validBlocksString );
  }

  LesHouchesAccordBlockEntryManager::LesHouchesAccordBlockEntryManager(
                                 std::set< std::string > const& validBlocksSet,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                         std::string const& fixedScaleArgument,
                                           std::string const& maximumScaleType,
                                    std::string const& maximumScaleArgument ) :
    LagrangianParameterManager(),
    numberOfDistinctActiveParameters( 0 ),
    activeParametersToIndices(),
    referenceSafeActiveParameters(),
    referenceUnsafeActiveParameters(),
    validBlocks( validBlocksSet ),
    lhaParser(),
    minimumScaleType( minimumScaleType ),
    minimumScaleArgument( minimumScaleArgument ),
    fixedScaleType( fixedScaleType ),
    fixedScaleArgument( fixedScaleArgument ),
    maximumScaleType( maximumScaleType ),
    maximumScaleArgument( maximumScaleArgument )
  {
    // This constructor is just an initialization list.
  }

  LesHouchesAccordBlockEntryManager::LesHouchesAccordBlockEntryManager(
                                             std::string const& xmlFileName ) :
    LagrangianParameterManager(),
    numberOfDistinctActiveParameters( 0 ),
    activeParametersToIndices(),
    referenceSafeActiveParameters(),
    referenceUnsafeActiveParameters(),
    validBlocks(),
    lhaParser(),
    minimumScaleType( "FixedNumber" ),
    minimumScaleArgument( "1.0" ),
    fixedScaleType( "FixedNumber" ),
    fixedScaleArgument( "1.0" ),
    maximumScaleType( "FixedNumber" ),
    maximumScaleArgument( "1.0" )
  {
    LHPC::RestrictedXmlParser xmlParser;
    xmlParser.OpenRootElementOfFile( xmlFileName );
    std::string renormalizationScaleChoices;
    while( xmlParser.ReadNextElement() )
    {
      if( xmlParser.CurrentName() == "RenormalizationScaleChoices" )
      {
        renormalizationScaleChoices = xmlParser.TrimmedCurrentBody();
      }
      else if( xmlParser.CurrentName() == "ValidBlocks" )
      {
        ParseValidBlocks( xmlParser.TrimmedCurrentBody() );
      }
    }
    xmlParser.CloseFile();
    xmlParser.LoadString( renormalizationScaleChoices );
    while( xmlParser.ReadNextElement() )
    {
      if( xmlParser.CurrentName() == "MinimumScaleBound" )
      {
        ParseScaleTypeAndArgument( xmlParser.TrimmedCurrentBody(),
                                   minimumScaleType,
                                   minimumScaleArgument );
      }
      else if( xmlParser.CurrentName() == "FixedScaleChoice" )
      {
        ParseScaleTypeAndArgument( xmlParser.TrimmedCurrentBody(),
                                   fixedScaleType,
                                   fixedScaleArgument );
      }
      else if( xmlParser.CurrentName() == "MaximumScaleBound" )
      {
        ParseScaleTypeAndArgument( xmlParser.TrimmedCurrentBody(),
                                   maximumScaleType,
                                   maximumScaleArgument );
      }
    }
  }

  LesHouchesAccordBlockEntryManager::~LesHouchesAccordBlockEntryManager()
  {
    for( size_t deletionIndex( 0 );
         deletionIndex < referenceSafeActiveParameters.size();
         ++deletionIndex )
    {
      delete referenceSafeActiveParameters[ deletionIndex ];
    }
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
    size_t const openBracket( trimmedVariable.find( '[' ) );
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
  std::string LesHouchesAccordBlockEntryManager::AsDebuggingString() const
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
    stringBuilder << "referenceSafeActiveParameters = { " << std::endl;
    for( std::vector< LhaBlockEntryInterpolator* >::const_iterator
         interpolatedParameter( referenceSafeActiveParameters.begin() );
         interpolatedParameter < referenceSafeActiveParameters.end();
         ++interpolatedParameter )
    {
      stringBuilder
      << (*interpolatedParameter)->AsDebuggingString() << std::endl;
    }
    stringBuilder
    << "}" << std::endl << "referenceUnsafeActiveParameters = { " << std::endl;
    for( std::vector< LhaBlockEntryInterpolator >::const_iterator
         interpolatedParameter( referenceUnsafeActiveParameters.begin() );
         interpolatedParameter < referenceUnsafeActiveParameters.end();
         ++interpolatedParameter )
    {
      stringBuilder
      << interpolatedParameter->AsDebuggingString() << std::endl;
    }
    stringBuilder << "}" << std::endl << "validBlocks = { ";
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
    std::vector< double > fixedScaleparameters;
    ParameterValues( log( AppropriateSingleFixedScale() ),
                     fixedScaleparameters );
    stringBuilder
    << "ParameterValues( log( AppropriateSingleFixedScale(), x ) makes x = { ";
    for( std::vector< double >::const_iterator
         fixedScaleParameter( fixedScaleparameters.begin() );
         fixedScaleParameter != fixedScaleparameters.end();
         ++fixedScaleParameter )
    {
      stringBuilder << *fixedScaleParameter << ", ";
    }
    stringBuilder << "}" << std::endl;
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
