/*
 * ParsingUtilities.hpp
 *
 *  Created on: Nov 16, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LHPC_PARSINGUTILITIES_HPP_
#define LHPC_PARSINGUTILITIES_HPP_

#include <string>
#include <sstream>

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
    // If there are any more characters in validBlocks which are not in
    // blockSeparators, we have at least one substring to add.
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

}

#endif /* LHPC_PARSINGUTILITIES_HPP_ */
