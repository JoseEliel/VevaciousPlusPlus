/*
 * RunningParameterManager.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  RunningParameterManager::RunningParameterManager() :
    parameterFunctionoidPointers(),
    parameterFunctionoidMap(),
    functionoidFinder(),
    slhaParser(),
    slhaBlockPointers(),
    slhaBlockMap(),
    slhaBlockFinder(),
    slhaAliasMap(),
    slhaAliasFinder()
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

    for( int deletionIndex( parameterFunctionoidPointers.size() - 1 );
         deletionIndex >= 0;
         --deletionIndex )
    {
      delete parameterFunctionoidPointers.at( deletionIndex );
    }
    for( int deletionIndex( slhaBlockPointers.size() - 1 );
         deletionIndex >= 0;
         --deletionIndex )
    {
      delete slhaBlockPointers.at( deletionIndex );
    }
  }


  // This adds to the list of SLHA blocks that will be read in when updating
  // running parameters. An alias for the block can be provided.
  void
  RunningParameterManager::AddValidSlhaBlock( std::string slhaBlock,
                                              std::string const& aliasString )
  {
    BOL::StringParser::transformToUppercase( slhaBlock );
    bool notAlreadyDeclared( slhaBlockMap.insert(
                     std::pair< std::string, RestrictedSlhaBlock* >( slhaBlock,
                                                             NULL ) ).second );
    if( !notAlreadyDeclared )
    {
      throw std::runtime_error(
              "Multiple declarations of the same SLHA block in model file!" );
    }
    if( !(aliasString.empty()) )
    {
      notAlreadyDeclared = slhaAliasMap.insert(
                            std::pair< std::string, std::string >( aliasString,
                                                        slhaBlock ) ).second;
      if( !notAlreadyDeclared )
      {
        throw std::runtime_error(
                        "Same alias for multiple SLHA blocks in model file!" );
      }
    }
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

    size_t wordEnd( parameterString.find_first_of( '=' ) );
    if( ( wordEnd == std::string::npos )
        ||
        ( parameterString[ parameterString.size() - 1 ] != ']' ) )
    {
      throw
      std::runtime_error( "Model file had malformed derived parameter!" );
    }
    std::string aliasString( BOL::StringParser::trimFromFrontAndBack(
                                                     parameterString.substr( 0,
                                                                 wordEnd ) ) );
    size_t wordStart( wordEnd );
    wordEnd = parameterString.find_first_of( '[',
                                             wordStart );
    // The extra (+/-) 1 and 2 in the substr arguments just remove the '[' and
    // ']' from being included.
    std::string typeString( BOL::StringParser::trimFromFrontAndBack(
                                     parameterString.substr( ( wordStart + 1 ),
                                             ( wordEnd - wordStart - 1 ) ) ) );
    std::string bracketedString( parameterString.substr( ( wordEnd + 1 ),
                                 ( parameterString.size() - wordEnd - 2 ) ) );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "aliasString = \"" << aliasString << "\", typeString = \"" << typeString
    << "\", bracketedString = \"" << bracketedString << "\"";
    std::cout << std::endl;/**/

    ParameterFunctionoid* createdFunctionoid( NULL );
    bool notAlreadyDeclared( true );
    if( typeString.compare( "SLHAELEMENT" ) == 0 )
    {
      wordEnd = bracketedString.find_first_of( ',' );
      createdFunctionoid
      = GetSlhaFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                     bracketedString.substr( 0,
                                                                   wordEnd ) ),
                    FormatIndexBracketContent( bracketedString.substr( wordEnd ) ) );
      if( createdFunctionoid == NULL )
      {
        throw std::runtime_error( "Derived SLHAELEMENT parameter referenced"
                                  " undeclared SLHA block in model file!" );
      }
    }
    notAlreadyDeclared = parameterFunctionoidMap.insert(
                  std::pair< std::string, ParameterFunctionoid* >( aliasString,
                                                 createdFunctionoid ) ).second;
    if( !notAlreadyDeclared )
    {
      throw std::runtime_error(
                 "Same alias for multiple derived parameters in model file!" );
    }
  }

  // This returns a pointer to the ParameterFunctionoid for the block named
  // blockName with indices given by indexString, making a new one if there
  // is none already. If the block name was invalid though, nothing is made
  // and NULL is returned.
  ParameterFunctionoid*
  RunningParameterManager::GetSlhaFunctionoid( std::string const& blockName,
                                               std::string const& indexString )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "RunningParameterManager::GetSlhaFunctionoid( \"" << blockName
    << "\", \"" << indexString << "\" ) called.";
    std::cout << std::endl;/**/

    ParameterFunctionoid* returnPointer( NULL );

    // First we check to see if blockName was an alias:
    std::string const* properBlockName( &blockName );
    slhaAliasFinder = slhaAliasMap.find( blockName );
    if( slhaAliasFinder != slhaAliasMap.end() )
    {
      properBlockName = &(slhaAliasFinder->second);
    }


    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "*properBlockName = " << *properBlockName;
    std::cout << std::endl;/**/

    slhaBlockFinder = slhaBlockMap.find( *properBlockName );
    if( slhaBlockFinder != slhaBlockMap.end() )
    {
      // If *properBlockName was a valid declared SLHA block, we form the
      // string which the functionoid should be mapped to.
      std::string lookupString( *properBlockName );
      lookupString.append( "[" );
      lookupString.append( indexString );
      lookupString.append( "]" );
      // We now check to see if the proper mapping string already exists
      // (which may happen if the original blockName was an alias).
      functionoidFinder = parameterFunctionoidMap.find( lookupString );
      if( functionoidFinder != parameterFunctionoidMap.end() )
      {
        // If the mapping did exist for the alias, we insert the alias into the
        // map before returning the pointer.
        lookupString = blockName;
        lookupString.append( "[" );
        lookupString.append( indexString );
        lookupString.append( "]" );
        parameterFunctionoidMap.insert(
                 std::pair< std::string, ParameterFunctionoid* >( lookupString,
                                                 functionoidFinder->second ) );
        return functionoidFinder->second;
      }

      // Otherwise we need to make a new functionoid. We make the functionoid
      // before its SLHA block because it sets up the number of indices which
      // is needed for the block constructor.
      SlhaFunctionoid* slhaFunctionoid = new SlhaFunctionoid( indexString );
      parameterFunctionoidPointers.push_back( slhaFunctionoid );
      parameterFunctionoidMap.insert(
                 std::pair< std::string, ParameterFunctionoid* >( lookupString,
                                                           slhaFunctionoid ) );
      slhaFunctionoidPointers.push_back( slhaFunctionoid );
      if( slhaBlockFinder->second == NULL )
      {
        slhaBlockFinder->second = new RestrictedSlhaBlock( *properBlockName,
                                            slhaFunctionoid->NumberOfIndices(),
                                                           0.0 );
        slhaParser.registerBlock( *(slhaBlockFinder->second) );
        slhaBlockPointers.push_back( slhaBlockFinder->second );
      }
      slhaFunctionoid->SetBlockPointer( slhaBlockFinder->second );
      returnPointer = slhaFunctionoid;
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "slhaBlockMap:";
    for( std::map< std::string, RestrictedSlhaBlock* >::iterator
         whichPair( slhaBlockMap.begin() );
         whichPair != slhaBlockMap.end();
         ++whichPair )
    {
      std::cout << "\"" << whichPair->first << "\" -> " << whichPair->second;
      if( whichPair->second != NULL )
      {
        std::cout << " (\"" << whichPair->second->getName() << "\")";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl << "parameterFunctionoidMap:" << std::endl;
    for( std::map< std::string, ParameterFunctionoid* >::iterator
         whichPair( parameterFunctionoidMap.begin() );
         whichPair != parameterFunctionoidMap.end();
         ++whichPair )
    {
      std::cout << "\"" << whichPair->first << "\" -> " << whichPair->second;
      if( whichPair->second != NULL )
      {
        std::cout << " (\"" << whichPair->second->AsString() << "\")";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;/**/


    return returnPointer;
  }

} /* namespace VevaciousPlusPlus */
