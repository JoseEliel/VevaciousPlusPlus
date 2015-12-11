/*
 * ParsingUtilities.hpp
 *
 *  Created on: Nov 16, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef LHPC_PARSINGUTILITIES_HPP_
#define LHPC_PARSINGUTILITIES_HPP_

#include <string>
#include <sstream>
#include <cstdlib>

namespace LHPC
{
  // This class is for parsing strings. The boost libraries probably do it all
  // and better, but LHPC is not meant to depend on boost.
  class ParsingUtilities
  {
  public:
    static std::string WhitespaceChars() { return " \t"; }

    static std::string NewlineChars() { return "\n\r"; }

    static std::string WhitespaceAndNewlineChars() { return " \t\n\r"; }

    static std::string LowercaseAlphabetChars()
    { return "abcdefghijklmnopqrstuvwxyz"; }

    static std::string UppercaseAlphabetChars()
    { return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"; }

    static std::string DigitChars() { return "0123456789"; }

    static void
    TransformToUppercase( std::string& stringToTransform );

    static std::string TrimFromFrontAndBack( std::string const& stringToTrim,
                                             std::string const& charsToTrim );

    static std::string TrimFromFront( std::string const& stringToTrim,
                                      std::string const& charsToTrim )
    { return
      stringToTrim.substr( stringToTrim.find_first_not_of( charsToTrim ),
                           std::string::npos ); }

    static std::string TrimFromBack( std::string const& stringToTrim,
                                     std::string const& charsToTrim )
    { return stringToTrim.substr( 0,
                              stringToTrim.find_last_not_of( charsToTrim ) ); }

    static std::string
    TrimWhitespaceFromFrontAndBack( std::string const& stringToTrim )
    { return TrimFromFrontAndBack( stringToTrim,
                                   WhitespaceAndNewlineChars() ); }

    // This splits stringToSplit into a vector of non-empty strings, which are
    // the substrings of all characters in stringToChange which do not appear
    // in separatorCharacters.
    static std::vector< std::string >
    SplitBySubstrings( std::string const& stringToChange,
                       std::string const& separatorCharacters );

    // This is just a shorthand for a commonly-used conversion.
    static long int BaseTenStringToInt( std::string const& baseTenString )
    { return strtol( baseTenString.c_str(),
                     NULL,
                     10 ); }

    // This is just a shorthand for a commonly-used conversion.
    static double StringToDouble( std::string const& baseTenString )
    { return strtod( baseTenString.c_str(),
                     NULL ); }

    // This parses indicesString as a set of integers separated by non-digit
    // characters and returns the set as a vector.
    static std::vector< int > ParseIndices( std::string const& indicesString );

    // This returns the position in contentLine of the first non-whitespace
    // character after the substring of contentLine which matches the indices
    // in indicesVector. If the indices are not matched, std::string::npos is
    // returned. (Also, if the line is nothing but the indices,
    // std::string::npos will be returned in this case as well.)
    static size_t StartOfMatchedContent( std::string const& contentLine,
                                     std::vector< int > const& indicesVector );

    // This returns the position in contentLine of the first non-whitespace
    // character after the substring of contentLine starting from startPosition
    // which matches the index in indexValue. If the index is not matched,
    // std::string::npos is returned. (Also, if the substring is nothing but
    // the index, std::string::npos will be returned in this case as well.) It
    // is assumed that all leading whitespace has already been skipped by
    // startPosition.
    static size_t MatchesIndex( std::string const& contentLine,
                                size_t startPosition,
                                int const indexValue );

    // This replaces all instances of characterToRemove in stringToModify with
    // characterToInsert.
    static void ReplaceAllCharacters( std::string& stringToModify,
                                      char const characterToRemove,
                                      char const characterToInsert );

    // This returns true if stringToSearch contains at least one instance of
    // characterToSeek.
    static bool CharacterIsInString( char const characterToSeek,
                                     std::string const& stringToSearch );

    // This resets the flags of streamToReset and sets its string to be
    // newString.
    static void ResetStringstream( std::stringstream& streamToReset,
                                   std::string const& newString );

    // This resets the flags of streamToReset and sets its string to be "".
    static void ResetStringstream( std::stringstream& streamToReset )
    { ResetStringstream( streamToReset,
                         "" ); }

    // This resets the flags of streamToReset and sets its string to be
    // newString.
    static void ResetStringstream( std::istringstream& streamToReset,
                                   std::string const& newString );

    // This returns the given double in the form "(1.234567 * 10^(-8))".
    static std::string
    FormatNumberForMathematica( double const numberToFormat );
  };





  inline void
  ParsingUtilities::TransformToUppercase( std::string& stringToTransform )
  {
    for( std::string::iterator stringCharacter( stringToTransform.begin() );
         stringCharacter != stringToTransform.end();
         ++stringCharacter )
    {
      *stringCharacter = std::toupper( *stringCharacter );
    }
  }

  inline std::string
  ParsingUtilities::TrimFromFrontAndBack( std::string const& stringToTrim,
                                          std::string const& charsToTrim )
  {
    size_t startPosition( stringToTrim.find_first_not_of( charsToTrim ) );
    if( startPosition == std::string::npos )
    {
      return "";
    }
    return stringToTrim.substr( startPosition,
                                ( stringToTrim.find_last_not_of( charsToTrim )
                                  + 1
                                  - startPosition ) );
  }

  // This splits stringToSplit into a vector of non-empty strings, which are
  // the substrings of all characters in stringToChange which do not appear
  // in separationCharacters.
  inline std::vector< std::string >
  ParsingUtilities::SplitBySubstrings( std::string const& stringToSplit,
                                      std::string const& separationCharacters )
  {
    std::vector< std::string > returnVector;
    size_t
    wordStart( stringToSplit.find_first_not_of( separationCharacters ) );
    size_t wordEnd( 0 );
    // If there are any more characters in stringToSplit which are not in
    // separationCharacters, we have at least one substring to add.
    while( wordStart != std::string::npos )
    {
      wordEnd = stringToSplit.find_first_of( separationCharacters,
                                             wordStart );
      returnVector.push_back( stringToSplit.substr( wordStart,
                                                   ( wordEnd - wordStart ) ) );
      if( wordEnd == std::string::npos )
      {
        break;
      }
      wordStart = stringToSplit.find_first_not_of( separationCharacters,
                                                   wordEnd );
    }
    return returnVector;
  }

  // This parses indicesString as a set of integers separated by non-digit
  // characters and returns the set as a vector.
  inline std::vector< int >
  ParsingUtilities::ParseIndices( std::string const& indicesString )
  {
    std::vector< int > indicesVector;
    size_t wordStart( indicesString.find_first_of(
                                            ParsingUtilities::DigitChars() ) );
    size_t wordEnd( 0 );
    while( wordStart != std::string::npos )
    {
      wordEnd
      = indicesString.find_first_not_of( ParsingUtilities::DigitChars(),
                                         wordStart );
      indicesVector.push_back( BaseTenStringToInt( indicesString.substr(
                                                                     wordStart,
                                                 ( wordEnd - wordStart ) ) ) );
      if( wordEnd == std::string::npos )
      {
        break;
      }
      wordStart = indicesString.find_first_of( ParsingUtilities::DigitChars(),
                                               wordEnd );
    }
    return indicesVector;
  }

  // This returns the position in contentLine of the first non-whitespace
  // character after the substring of contentLine which matches the indices
  // in indicesVector. If the indices are not matched, std::string::npos is
  // returned. (Also, if the line is nothing but the indices,
  // std::string::npos will be returned in this case as well.)
  inline size_t
  ParsingUtilities::StartOfMatchedContent( std::string const& contentLine,
                                      std::vector< int > const& indicesVector )
  {
    size_t contentStart( contentLine.find_first_not_of(
                                     ParsingUtilities::WhitespaceChars() ) );
    for( std::vector< int >::const_iterator
         indexValue( indicesVector.begin() );
         indexValue != indicesVector.end();
         ++indexValue )
    {
      contentStart = MatchesIndex( contentLine,
                                   contentStart,
                                   *indexValue );
      if( contentStart == std::string::npos )
      {
        return std::string::npos;
      }
    }
    return contentStart;
  }

  // This returns the position in contentLine of the first non-whitespace
  // character after the substring of contentLine starting from startPosition
  // which matches the index in indexValue. If the index is not matched,
  // std::string::npos is returned. (Also, if the substring is nothing but
  // the index, std::string::npos will be returned in this case as well.) It
  // is assumed that all leading whitespace has already been skipped by
  // startPosition.
  inline size_t ParsingUtilities::MatchesIndex( std::string const& contentLine,
                                                size_t startPosition,
                                                int const indexValue )
  {
    size_t
    indexEnd( contentLine.find_first_of( ParsingUtilities::WhitespaceChars(),
                                         startPosition ) );
    if( indexValue == BaseTenStringToInt( contentLine.substr( startPosition,
                                           ( indexEnd - startPosition ) ) ) )
    {
      return
      contentLine.find_first_not_of( ParsingUtilities::WhitespaceChars(),
                                     indexEnd );
    }
    else
    {
      return std::string::npos;
    }
  }

  // This replaces all instances of characterToRemove in stringToModify with
  // characterToInsert.
  inline void
  ParsingUtilities::ReplaceAllCharacters( std::string& stringToModify,
                                          char const characterToRemove,
                                          char const characterToInsert )
  {
    for( size_t characterIndex( 0 );
         characterIndex < stringToModify.size();
         ++characterIndex )
    {
      if( stringToModify[ characterIndex ] == characterToRemove )
      {
        stringToModify[ characterIndex ] = characterToInsert;
      }
    }
  }

  // This returns true if stringToSearch contains at least one instance of
  // characterToSeek.
  inline bool
  ParsingUtilities::CharacterIsInString( char const characterToSeek,
                                         std::string const& stringToSearch )
  {
    for( std::string::const_iterator stringCharacter(stringToSearch.begin());
         stringCharacter != stringToSearch.end();
         ++stringCharacter )
    {
      if( *stringCharacter == characterToSeek )
      {
        return true;
      }
    }
    return false;
  }

  // This resets the flags of streamToReset and sets its string to be
  // newString.
  inline void
  ParsingUtilities::ResetStringstream( std::stringstream& streamToReset,
                                       std::string const& newString )
  {
    streamToReset.clear();
    streamToReset.str( newString );
  }

  // This resets the flags of streamToReset and sets its string to be
  // newString.
  inline void
  ParsingUtilities::ResetStringstream( std::istringstream& streamToReset,
                                       std::string const& newString )
  {
    streamToReset.clear();
    streamToReset.str( newString );
  }

  // This returns the given double in the form "(1.234567 * 10^(-8))".
  inline std::string
  ParsingUtilities::FormatNumberForMathematica( double const numberToFormat )
  {
    std::stringstream stringBuilder;
    stringBuilder << '(' << numberToFormat << ')';
    std::string returnString( stringBuilder.str() );
    size_t exponentPosition( returnString.find_first_of( "eE" ) );
    if( exponentPosition == std::string::npos )
    {
      return returnString;
    }
    return ( returnString.replace( exponentPosition,
                                   1,
                                   "* 10^(" ) + ")" );
  }

} /* namespace LHPC */

#endif /* LHPC_PARSINGUTILITIES_HPP_ */
