/*
 * RestrictedXmlParser.hpp
 *
 *  Created on: Dec 8, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef RESTRICTEDXMLPARSER_HPP_
#define RESTRICTEDXMLPARSER_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include "ParsingUtilities.hpp"

namespace LHPC
{
  // Objects of this class parse out blocks of ASCII text from a string between
  // XML opening and closing tags, and return the text between the tags
  // (without interpreting it further as XML).
  // This class does not parse general XML text, and is not compliant with all
  // requirements for an XML parser.
  // Reference http://www.w3.org/TR/2008/REC-xml-20081126/ for XML
  // definition.
  //
  // This class is not properly compliant with XML in the following senses:
  // 1) It does not parse any prolog. If it is given an XML file to open, it
  //    looks for the first possible element, which is taken to be the first
  //    occurrence of '<' followed by any of the valid ASCII start characters
  //    for names (":" | [A-Z] | "_" | [a-z]). If comments or CDATA or ENTITY
  //    structures include such a substring, then the parsing will fail badly.
  //    Any text encountered before the assumed root element is stored
  //    literally as the prolog.
  //    If the class is passed a string rather than told to open a file, it
  //    assumes that the content is the content of an XML element and thus
  //    can contain only child elements, comments, CDATA, and miscellaneous
  //    text.
  // 2) It does not validate much. It validates some stuff, but not everything.
  // 3) It does not accept UTF-16.
  //
  // This class is compatible with UTF-8 in a sense, as the bytes will be read
  // correctly, but they will be stored in string objects, which are
  // effectively arrays of char data, so (depending on architecture in a way
  // that is beyond my expertise) will be stored as a series of bytes.
  // Conversion back to UTF-8 is beyond the scope of this class. The Boost
  // library provides string UTF conversion tools.
  // (http://www.boost.org/doc/libs/1_49_0/libs/locale/doc/html/charset_handling.html)
  // It's worth checking this out:
  // http://www.joelonsoftware.com/articles/Unicode.html on character encoding.
  //
  // Notes to self:
  // #x9 is '\t'
  // #xA is '\n'
  // #xD is '\r'
  // #x20 is ' '
  class RestrictedXmlParser
  {
  public:
    RestrictedXmlParser() : fileProlog( "" ),
                            rootName( "" ),
                            rootAttributes(),
                            currentName( "" ),
                            currentAttributes(),
                            currentBody( "" ),
                            xmlFileStream(),
                            xmlStringStream(),
                            xmlStream( NULL ),
                            currentChar( '?' ),
                            currentElementIsEmpty( true ) {}
    ~RestrictedXmlParser() { CloseFile(); }


    // This loads the given string as the content to parse into XML elements.
    void LoadString( std::string const& stringToParse )
    {
      ResetContent();
      ParsingUtilities::ResetStringstream( xmlStringStream,
                                           stringToParse );
      xmlStream = &xmlStringStream;
    }

    // This opens the file with name fileName, and reads until the first
    // element start tag is found, storing the prolog, and reads in that first
    // start tag, then stops. The prolog is accessible through FileProlog(),
    // and the content of the start tag is accessible through RootName() and
    // RootAttributes(). The rest of the root element is treated as the
    // content to parse into XML elements, as stringToParse would be in
    // LoadString(...) above, but does not load the whole file into memory,
    // rather just streams the chars from the fstream. Note that the file will
    // remain open after calling this function, until CloseFile() is called.
    // (The destructor calls CloseFile() in case the user forgot, CloseFile()
    // is also called at the start of this function to ensure that the last
    // file (if any was open) really does get shut.)
    void OpenRootElementOfFile( std::string const& fileName )
    {
      CloseFile();
      xmlFileStream.open( fileName.c_str() );
      xmlStream = &xmlFileStream;
      ReadPrologAndOpenRootElement();
    }

    // This closes the file which was being parsed, if there is one open.
    void CloseFile()
    {
      if( xmlFileStream.is_open() )
      {
        xmlFileStream.close();
      }
      xmlFileStream.clear();
      xmlStream = &xmlStringStream;
    }

    // This opens the root element of the file as done by
    // OpenRootElementOfFile(...), then reads in the entire content of the
    // root element, with the consequence that CurrentName() and
    // CurrentAttributes() are the same as RootName() and RootAttributes()
    // respectively, and CurrentBody() is the full content of the root element.
    // It does close the file as part of the function, so CloseFile() is not
    // necessary afterwards.
    void ReadAllOfRootElementOfFile( std::string const& fileName )
    {
      OpenRootElementOfFile( fileName );
      RecordToEndOfElement( rootName );
      CloseFile();
    }

    // This moves forward in the current content to parse until the next
    // element start tag is found, then that element is read in until its end
    // tag is found, and the character just after that is where the next call
    // of this function will resume. The element just parsed is accessible
    // through CurrentName() for its name, CurrentAttributes() for its
    // attributes, and CurrentBody() for the text between the start and end
    // tags. The return value is true if an element was read in, or false if no
    // more valid elements could be found. (If no further elements are
    // found, the next call of this function will do nothing, as the next
    // character will just be the EOF character.)
    bool ReadNextElement();// { NEEDS DOING! }

    // This resets the parser so that ReadNextElement() will start from the
    // beginning of the content to be parsed (which is either the string given
    // to LoadString(...) or the body of the root element of the file given to
    // OpenRootElementOfFile(...); if called after calling
    // ReadAllOfRootElementOfFile(...), nothing happens).
    void ReturnToBeginningOfText()
    {
      ResetContent();
      xmlStream->clear();
      xmlStream->seekg( std::ios::beg );
    }

    // This returns the name of the element just read by the last call of
    // ReadNextElement().
    std::string const& CurrentName() { return currentName; }

    // This returns the attributes of the element just read by the last call of
    // ReadNextElement().
    std::map< std::string, std::string > const& CurrentAttributes()
    { return currentAttributes; }

    // This returns the body of the element just read by the last call of
    // ReadNextElement().
    std::string const& CurrentBody() { return currentBody; }

    // This returns the text of the file opened by OpenRootElementOfFile(...)
    // or ReadAllOfRootElementOfFile(...) up to the beginning of the root
    // element.
    std::string const& FileProlog() { return fileProlog; }

    // This returns the name of the root element of the file opened by
    // OpenRootElementOfFile(...) or ReadAllOfRootElementOfFile(...).
    std::string const& RootName() { return rootName; }

    // This returns the attributes of the root element of the file opened by
    // OpenRootElementOfFile(...) or ReadAllOfRootElementOfFile(...).
    std::map< std::string, std::string > const& RootAttributes()
    { return rootAttributes; }


  protected:
    typedef std::map< std::string, std::string > AttributeMap;
    static std::string const AllowedWhitespaceChars() { return " \t\r\n"; }
    static std::string const AllowedQuoteChars() { return "\'\""; }
    static std::string const AllowedNameStartChars()
    { return ( ParsingUtilities::UppercaseAlphabetChars()
               + ParsingUtilities::LowercaseAlphabetChars()
               + ":_" ); }
    static std::pair< std::string, std::string > const
    CommentDelimiter() { return std::make_pair( "<!--", "-->" ); }
    static std::pair< std::string, std::string > const
    PiDelimiter() { return std::make_pair( "<?", "?>" ); }
    static std::pair< std::string, std::string > const
    DoctypeDelimiter() { return std::make_pair( "<!DOCTYPE", ">" ); }
    static std::pair< std::string, std::string > const
    CdataDelimiter() { return std::make_pair( "<![CDATA[", "]]>" ); }

    // Data members stored from parsing:
    std::string fileProlog;
    std::string rootName;
    AttributeMap rootAttributes;
    std::string currentName;
    AttributeMap currentAttributes;
    std::string currentBody;

    // Data members used by parsing:
    std::ifstream xmlFileStream;
    std::istringstream xmlStringStream;
    std::istream* xmlStream;
    char currentChar;
    bool currentElementIsEmpty;

    // This sets the various recording data to the values they should have
    // before reading in some text.
    void ResetContent()
    {
      fileProlog.assign( "" );
      rootName.assign( "" );
      rootAttributes.clear();
      currentName.assign( "" );
      currentAttributes.clear();
      currentBody.assign( "" );
    }

    // This resets the content data members, then reads in everything up to
    // (but not including) the '<' of the first valid opening tag of an
    // element, then reads the name in that tag into rootName, and the
    // attributes into rootAttributes, and leaves xmlStream ready to read in
    // the first character after the '>' which closes the opening tag of the
    // root element.
    void ReadPrologAndOpenRootElement()
    {
      ResetContent();
      std::stringstream prologStream;
      std::stringstream tagStream;
      if( ReadToNextTagOpener( &prologStream ) )
      {
        ReadStartTag( rootName,
                      &rootAttributes,
                      tagStream );
      }
      throw std::runtime_error( "No root element found in file!" );
    }

    // This puts characters from xmlStream into nonMarkupDestination if not
    // NULL, or discards the characters if nonMarkupDestination is NULL, until
    // it reads in the first instance of '<' followed by an allowed element
    // name start character. Neither the '<' nor the starting character are put
    // into nonMarkupDestination, but the starting character is left in
    // currentChar. True is returned if such a sequence is found, false if not
    // and no more characters can be read.
    bool ReadToNextTagOpener( std::istream* nonMarkupDestination )
    {
      while( ReadToNextHalt( '<',
                             nonMarkupDestination,
                             NULL )
             &&
             xmlStream->get( currentChar ).good() )
      {
        // We break from the loop with true as soon as we find a '<' followed
        // by an allowed name start character.
        if( ParsingUtilities::CharacterIsInString( currentChar,
                                                   AllowedNameStartChars() ) )
        {
          return true;
        }
        else if( nonMarkupDestination != NULL )
        {
          (*nonMarkupDestination) << '<' << currentChar;
        }
      }
      // If we did not break from the loop, then a tag opener was not found and
      // xmlStream is finished.
      return false;
    }

    // This puts characters from xmlStream into destinationWithoutHalt if it is
    // not NULL and destinationIncludingHalt if it is not NULL (if both are
    // NULL, the characters are discarded) until it reads in the first instance
    // of haltCharacter (which is put into destinationIncludingHalt but not
    // into destinationWithoutHalt).
    bool ReadToNextHalt( char const haltCharacter,
                         std::istream* destinationWithoutHalt,
                         std::istream* destinationIncludingHalt )
    {
      while( ReadCharacter( destinationIncludingHalt )
             &&
             ( currentChar != haltCharacter ) )
      {
        if( destinationWithoutHalt != NULL )
        {
          (*destinationWithoutHalt) << currentChar;
        }
      }
      return ( currentChar == haltCharacter );
    }

    // This puts the next character from xmlStream into currentChar, then puts
    // it into destinationForReadCharacters if it is not NULL, then returns
    // true, unless no character could be read, in which case false is
    // returned.
    bool ReadCharacter( std::istream* destinationForReadCharacters )
    {
      if( xmlStream->get( currentChar ).good() )
      {
        if( destinationForReadCharacters != NULL )
        {
          (*destinationForReadCharacters) << currentChar;
        }
        return true;
      }
      else
      {
        return false;
      }
    }

    // This tries to read a start tag, assuming that '<' was read before
    // currentChar, and that currentChar is the first character of a valid
    // element name. It sets currentElementIsEmpty to be true if the tag was
    // actually an empty-element tag rather than a start tag. All the
    // characters of the tag are passed into tagRecord. If it fails to read a
    // valid start tag, it throws an exception.
    void ReadStartTag( std::string& nameDestination,
                            AttributeMap* attributeDestination,
                            std::stringstream& tagRecord )
    {
      ParsingUtilities::ResetStringstream( tagRecord );
      tagRecord << '<' << currentChar;
      std::stringstream bufferStream;
      bufferStream << currentChar;

      // We read to the first character that marks the end of the element name.
      if( !(ReadToNextHalt( ( AllowedWhitespaceChars() + ">/" ),
                            &bufferStream,
                            &tagRecord )) )
      {
        throw std::runtime_error(
               "Failed to find whitespace or end of tag after element name!" );
      }
      nameDestination = bufferStream.str();
      if( attributeDestination != NULL )
      {
        attributeDestination->clear();
      }
      CloseStartTag( attributeDestination,
                     tagRecord );
    }

    // This puts characters from xmlStream into destinationWithoutHalt if it is
    // not NULL and destinationIncludingHalt if it is not NULL (if both are
    // NULL, the characters are discarded) until it reads in the first instance
    // of any of the characters in haltCharacters (which is put into
    // destinationIncludingHalt but not into destinationWithoutHalt).
    bool ReadToNextHalt( std::string const& haltCharacters,
                         std::istream* destinationWithoutHalt,
                         std::istream* destinationIncludingHalt )
    {
      while( ReadCharacter( destinationIncludingHalt )
             &&
             !(ParsingUtilities::CharacterIsInString( currentChar,
                                                      haltCharacters )) )
      {
        if( destinationWithoutHalt != NULL )
        {
          (*destinationWithoutHalt) << currentChar;
        }
      }
      return ParsingUtilities::CharacterIsInString( currentChar,
                                                    haltCharacters );
    }

    // This reads xmlStream until the end of the start or empty-element tag is
    // reached, parsing attributes along the way into attributeDestination,
    // throwing an exception if it does not reach the end of the tag without
    // the stream ending or finding malformed XML. All the characters of the
    // tag get put into tagRecord if it is not NULL.
    void CloseStartTag( AttributeMap* attributeDestination,
                        std::stringstream& tagRecord )
    {
      if( !(SkipWhitespace( &tagRecord )) )
      {
        std::stringstream errorBuilder;
        errorBuilder << "";
        throw std::runtime_error( errorBuilder.str() );
      }
      // Now currentChar is the end of the start tag ('>'), or the first
      // character of the end of an empty element tag ('/' in "/>") (or the
      // tag is malformed with a '/' out of place), or is the first character
      // of an attribute name.
      if( currentChar == '>' )
      {
        currentElementIsEmpty = false;
      }
      else if( currentChar == '/' )
      {
        if( ReadCharacter( &tagRecord )
            &&
            ( currentChar == '>' ) )
        {
          currentElementIsEmpty = true;
        }
      }
      else
      {
        ParseAttribute( attributeDestination,
                        tagRecord );
        CloseStartTag( attributeDestination,
                       tagRecord );
      }
    }

    // This reads characters from xmlStream into currentChar until the first
    // non-whitespace character is read in, returning true unless the stream
    // ended without any non-whitespace character being read in. If
    // destinationForReadCharacters is not NULL, all read characters (including
    // the first non-whitespace character) are put into it.
    bool SkipWhitespace( std::istream* destinationForReadCharacters )
    {
      // If the name is followed by whitespace, we keep going to the first
      // non-whitespace character, discarding all the whitespace characters
      // along the way.
      while( ParsingUtilities::CharacterIsInString( currentChar,
                                                   AllowedWhitespaceChars() ) )
      {
        // If we run out of characters from the stream before finding a
        // non-whitespace character, we return false. Evaluating the
        // conditional does work.
        if( !(ReadCharacter( destinationForReadCharacters )) )
        {
          return false;
        }
      }
      return true;
    }

    // This parses an attribute from xmlStream assuming that the first
    // character of the attribute name is already in currentChar, and puts the
    // attribute into attributeDestination. It throws an exception if the XML
    // is malformed.
    void ParseAttribute( AttributeMap* attributeDestination,
                         std::stringstream& tagRecord )
    {
      // When this function is called, currentChar is the first character of
      // the name of an attribute.
      std::stringstream nameStream;
      nameStream << currentChar;
      std::stringstream valueStream;
      if( !( ReadToNextHalt( ( AllowedWhitespaceChars() + "=" ),
                             &nameStream,
                             &tagRecord )
             &&
             SkipWhitespace( &tagRecord )
             &&
             ( currentChar == '=' )
             &&
             ReadCharacter( &tagRecord )
             &&
             SkipWhitespace( &tagRecord )
             &&
             ( ( currentChar == '\'' ) || ( currentChar == '\"' ) )
             &&
             ReadToNextHalt( currentChar,
                             &nameStream,
                             &valueStream ) ) )
      {
        throw std::runtime_error( "Could not parse an attribute correctly!" );
      }
      (*attributeDestination)[ nameStream.str() ] = valueStream.str();
    }

    // This puts everything except comments from xmlStream into currentBody up
    // to the end tag for the element named elementName. If any nested start
    // tags are found for child elements with this name, the text is read in
    // until each child of that name is closed and then to the next end tag for
    // that name. It does not validate that any other tags or elements are
    // correctly nested and so on. It throws an exception if it cannot find the
    // end of the element.
    void RecordToEndOfElement( std::string const& elementName )
    {
      unsigned int numberOfUnclosedElementsOfGivenName( 1 );
      std::stringstream contentStream;
      std::stringstream tagStream;
      std::string innerName;
      while( ReadToNextHalt( '<',
                             &contentStream,
                             NULL ) )
      {
        if( !(xmlStream->get( currentChar ).good()) )
        {
          std::stringstream errorBuilder;
          errorBuilder
          << "Could not find end tag for element <" << elementName << ">";
          throw std::runtime_error( errorBuilder.str() );
        }

        if( currentChar == '?' )
        {
          CloseQuestionMark( &contentStream );
        }
        else if( currentChar == '!' )
        {
          CloseExclamationMark( &contentStream );
        }
        else if( ParsingUtilities::CharacterIsInString( currentChar,
                                                    AllowedNameStartChars() ) )
        {
          ReadStartTag( innerName,
                        NULL,
                        tagStream );
          contentStream << tagStream.str();
          if( !currentElementIsEmpty
              &&
              ( innerName == elementName ) )
          {
            ++numberOfUnclosedElementsOfGivenName;
          }
        }
        else if( currentChar == '/' )
        {
          ReadEndTag( innerName,
                      tagStream );
          if( innerName == elementName )
          {
            --numberOfUnclosedElementsOfGivenName;
            if( numberOfUnclosedElementsOfGivenName == 0 )
            {
              currentBody = contentStream.str();
              return;
            }
          }
          contentStream << tagStream.str();
        }
      }
    }

    // This tries to read an end tag, assuming that "</" was the last pair of
    // characters read in before this function was called, putting the name of
    // the element into nameDestination. The characters read in by this
    // function are passed to tagRecord. It throws an exception if it could not
    // read a valid end tag.
    void ReadEndTag( std::string& nameDestination,
                     std::stringstream& tagRecord )
    {
      tagRecord << "</";
      std::stringstream nameStream;
      if( !( ReadToNextHalt( ( AllowedWhitespaceChars() + ">" ),
                             &nameStream,
                             &tagRecord )
             &&
             SkipWhitespace( &tagRecord )
             &&
             ( currentChar == '>' ) ) )
      {
        throw std::runtime_error( "Could not close a valid end tag!" );
      }
      nameDestination = nameStream.str();
    }

    // This puts "<?" into destinationForReadCharacters if it is not NULL, then
    // reads in characters from xmlStream until the next "?>", placing
    // them all into destinationForReadCharacters if it is not NULL. It throws
    // an exception if it fails to read in "?>".
    void CloseQuestionMark( std::istream* destinationForReadCharacters )
    {
      if( destinationForReadCharacters != NULL )
      {
        destinationForReadCharacters << "<?";
      }
      while( ReadToNextHalt( '?',
                             NULL,
                             destinationForReadCharacters )
             &&
             ReadCharacter( destinationForReadCharacters ) )
      {
        if( currentChar == '>' )
        {
          return;
        }
      }
      throw std::runtime_error( "Failed to close processing instruction!" );
    }

    // This assumes that the characters just read from xmlStream were "<!", and
    // determines if this is the start of a comment or a CDATA structure or
    // just malformed XML. If it is a comment, it reads until the end of the
    // comment, discarding all the characters of the comment, and true is
    // returned. If it is a CDATA structure, characters are read in until the
    // end of the structure, and all the characters read are put into
    // destinationForReadCharacters if it is not NULL, and true is returned.
    // If the XML was malformed (including being unable to close a comment or
    // CDATA structure), false is returned.
    void CloseExclamationMark( std::istream* destinationForReadCharacters )
    {
      if( !(xmlStream->get( currentChar ).good()) )
      {
        throw
        std::runtime_error( "Failed to close structure following \"<!\"" );
      }
      if( currentChar == '-' )
      {
        ReadComment();
      }
      else if( currentChar == '[' )
      {
        ReadCdata( destinationForReadCharacters );
      }
    }

    // This tries to read in characters (discarding them) from xmlStream until
    // a valid comment has been read, assuming that the previous characters
    // were "<!-". It throws an exception if it fails to close the comment.
    void ReadComment()
    {
      while( xmlStream->get( currentChar ).good() )
      {
        if( ( currentChar == '-' )
            &&
            xmlStream->get( currentChar ).good()
            &&
            ( currentChar == '-' )
            &&
            xmlStream->get( currentChar ).good()
            &&
            ( currentChar == '>' ) )
        {
          return;
        }
      }
      throw std::runtime_error( "Failed to close comment!" );
    }

    // This puts "<!" into destinationForReadCharacters if it is not NULL, then
    // reads in characters from xmlStream until the next "]]>", placing
    // them all into destinationForReadCharacters if it is not NULL. It throws
    // an exception if it fails to read in "]]>".
    void ReadCdata( std::istream* destinationForReadCharacters )
    {
      if( destinationForReadCharacters != NULL )
      {
        destinationForReadCharacters << "<[";
      }
      while( ReadToNextHalt( ']',
                             NULL,
                             destinationForReadCharacters )
             &&
             ReadCharacter( destinationForReadCharacters )
             &&
             ( currentChar == ']' )
             &&
             ReadCharacter( destinationForReadCharacters ) )
      {
        if( currentChar == '>' )
        {
          return;
        }
      }
      throw std::runtime_error( "Failed to close CDATA structure!" );
    }



    // This feeds characterToAppend to destinationForReadCharacters if it is
    // not NULL.
    void AppendToStream( std::istream* destinationForReadCharacters,
                         char const characterToAppend )
    {
      if( destinationForReadCharacters != NULL )
      {
        (*destinationForReadCharacters) << characterToAppend;
      }
    }

    // This reads in characters from xmlStream and compares them to
    // stringToMatch. It reads the characters one at a time, and if a character
    // read in does not match that in the correct place in stringToMatch, the
    // matching characters so far plus this last-read character are put into
    // nonMatchDestination (if it is not NULL), and false is returned. If
    // xmlStream provides stringToMatch.size() characters in a row which match
    // stringToMatch exactly, then true is returned (without any characters
    // being put into nonMatchDestination).
    bool ReadMatch( std::string const& stringToMatch,
                    std::istream* nonMatchDestination )
    {
      for( size_t characterIndex( 0 );
           characterIndex < stringToMatch.size();
           ++characterIndex )
      {
        if( !( xmlStream->get( currentChar ).good()
               &&
               ( currentChar == stringToMatch[ characterIndex ] ) ) )
        {
          if( nonMatchDestination != NULL )
          {
            nonMatchDestination << stringToMatch.substr( 0,
                                                         characterIndex )
            << currentChar;
          }
          return false;
        }
      }
      return true;
    }






    //
    bool DiscardToNextTag();

    //
    bool ParseTagName( std::string& nameDestination );

    //
    bool ParseAttributes(
                  std::map< std::string, std::string >& attributeDestination );


    static char const markupOpener;
    static char const markupCloser;
    static char const tagCloser;
    static std::string const allowedXmlWhitespaceChars;
    static std::string const allowedXmlQuoteChars;
    static std::pair< std::string, std::string > const commentDelimiter;
    static std::pair< std::string, std::string > const piDelimiter;
    static std::pair< std::string, std::string > const doctypeDelimiter;
    static std::pair< std::string, std::string > const cdataDelimiter;
  };

}

#endif /* RESTRICTEDXMLPARSER_HPP_ */
