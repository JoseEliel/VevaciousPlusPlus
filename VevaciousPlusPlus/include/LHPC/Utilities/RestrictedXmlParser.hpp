/*
 * RestrictedXmlParser.hpp
 *
 *  Created on: Dec 8, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef RESTRICTEDXMLPARSER_HPP_
#define RESTRICTEDXMLPARSER_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <stdexcept>

#include "Utilities/ParsingUtilities.hpp"

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
                            currentCharacter( '?' ),
                            currentElementIsEmpty( true ) {}

    ~RestrictedXmlParser() { CloseFile(); }


    // This loads the given string as the content to parse into XML elements.
    void LoadString( std::string const& stringToParse );

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
    void OpenRootElementOfFile( std::string const& fileName );

    // This closes the file which was being parsed, if there is one open.
    void CloseFile();

    // This opens the root element of the file as done by
    // OpenRootElementOfFile(...), then reads in the entire content of the
    // root element, with the consequence that CurrentName() and
    // CurrentAttributes() are the same as RootName() and RootAttributes()
    // respectively, and CurrentBody() is the full content of the root element.
    // It does close the file as part of the function, so CloseFile() is not
    // necessary afterwards.
    void ReadAllOfRootElementOfFile( std::string const& fileName );

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
    bool ReadNextElement();

    // This resets the parser so that ReadNextElement() will start from the
    // beginning of the content to be parsed (which is either the string given
    // to LoadString(...) or the body of the root element of the file given to
    // OpenRootElementOfFile(...); if called after calling
    // ReadAllOfRootElementOfFile(...), nothing happens).
    void ReturnToBeginningOfText();

    // This returns the name of the element just read by the last call of
    // ReadNextElement().
    std::string const& CurrentName() const { return currentName; }

    // This returns the attributes of the element just read by the last call of
    // ReadNextElement().
    std::map< std::string, std::string > const& CurrentAttributes() const
    { return currentAttributes; }

    // This returns the body of the element just read by the last call of
    // ReadNextElement().
    std::string const& CurrentBody() const { return currentBody; }

    // This returns a string which is the current body with leading and
    // trailing whitespace and newline characters removed.
    std::string TrimmedCurrentBody() const
    { return ParsingUtilities::TrimWhitespaceFromFrontAndBack( currentBody ); }

    // This returns the text of the file opened by OpenRootElementOfFile(...)
    // or ReadAllOfRootElementOfFile(...) up to the beginning of the root
    // element.
    std::string const& FileProlog() const { return fileProlog; }

    // This returns the name of the root element of the file opened by
    // OpenRootElementOfFile(...) or ReadAllOfRootElementOfFile(...).
    std::string const& RootName() const { return rootName; }

    // This returns the attributes of the root element of the file opened by
    // OpenRootElementOfFile(...) or ReadAllOfRootElementOfFile(...).
    std::map< std::string, std::string > const& RootAttributes() const
    { return rootAttributes; }


  protected:
    typedef std::map< std::string, std::string > AttributeMap;

    static std::string const AllowedWhitespaceChars() { return " \t\r\n"; }
    static std::string const AllowedQuoteChars() { return "\'\""; }
    static std::string const AllowedNameStartChars()
    { return ( ParsingUtilities::UppercaseAlphabetChars()
               + ParsingUtilities::LowercaseAlphabetChars()
               + ":_" ); }


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
    char currentCharacter;
    bool currentElementIsEmpty;


    // This sets the various recording data to the values they should have
    // before reading in some text.
    void ResetContent();

    // This resets the content data members, then reads in everything up to
    // (but not including) the '<' of the first valid opening tag of an
    // element, then reads the name in that tag into rootName, and the
    // attributes into rootAttributes, and leaves xmlStream ready to read in
    // the first character after the '>' which closes the opening tag of the
    // root element.
    void ReadPrologAndOpenRootElement();

    // This puts characters from xmlStream into nonMarkupDestination if not
    // NULL, or discards the characters if nonMarkupDestination is NULL, until
    // it reads in the first instance of '<' followed by an allowed element
    // name start character. Neither the '<' nor the starting character are put
    // into nonMarkupDestination, but the starting character is left in
    // currentChar. True is returned if such a sequence is found, false if not
    // and no more characters can be read.
    bool ReadToNextTagOpener( std::ostream* nonMarkupDestination );

    // This either reads in a processing instruction or a comment or a CDATA
    // section, passing relevant characters to contentStream, or just puts '<'
    // followed by currentCharacter into contentStream, depending on what
    // currentCharacter is.
    void TryToCloseNonTagMarkup( std::ostream* destinationForReadCharacters );

    // This puts "<?" into destinationForReadCharacters if it is not NULL, then
    // reads in characters from xmlStream until the next "?>", placing
    // them all into destinationForReadCharacters if it is not NULL. It throws
    // an exception if it fails to read in "?>".
    void CloseQuestionMark( std::ostream* destinationForReadCharacters );

    // This puts characters from xmlStream into destinationWithoutHalt if it is
    // not NULL and destinationIncludingHalt if it is not NULL (if both are
    // NULL, the characters are discarded) until it reads in the first instance
    // of haltCharacter (which is put into destinationIncludingHalt but not
    // into destinationWithoutHalt).
    bool ReadToNextHalt( char const haltCharacter,
                         std::ostream* destinationWithoutHalt,
                         std::ostream* destinationIncludingHalt );

    // This puts the next character from xmlStream into currentChar, then puts
    // it into destinationForReadCharacters if it is not NULL, then returns
    // true, unless no character could be read, in which case false is
    // returned.
    bool ReadCharacter( std::ostream* destinationForReadCharacter );

    // This assumes that the characters just read from xmlStream were "<!", and
    // determines if this is the start of a comment or a CDATA structure or
    // just malformed XML. If it is a comment, it reads until the end of the
    // comment, discarding all the characters of the comment, and true is
    // returned. If it is a CDATA structure, characters are read in until the
    // end of the structure, and all the characters read are put into
    // destinationForReadCharacters if it is not NULL, and true is returned.
    // If the XML was malformed (including being unable to close a comment or
    // CDATA structure), false is returned.
    void CloseExclamationMark( std::ostream* destinationForReadCharacters );

    // This tries to read in characters (discarding them) from xmlStream until
    // a valid comment has been read, assuming that the previous characters
    // were "<!-". It throws an exception if it fails to close the comment.
    void ReadComment();

    // This puts "<!" into destinationForReadCharacters if it is not NULL, then
    // reads in characters from xmlStream until the next "]]>", placing
    // them all into destinationForReadCharacters if it is not NULL. It throws
    // an exception if it fails to read in "]]>".
    void ReadCdata( std::ostream* destinationForReadCharacters );

    // This tries to read a start tag, assuming that '<' was read before
    // currentChar, and that currentChar is the first character of a valid
    // element name. It sets currentElementIsEmpty to be true if the tag was
    // actually an empty-element tag rather than a start tag. All the
    // characters of the tag are passed into tagRecord if it is not NULL
    // (clearing the content of tagRecord beforehand). If it fails to read a
    // valid start tag, it throws an exception.
    void ReadStartTag( std::string& nameDestination,
                       AttributeMap* attributeDestination,
                       std::stringstream* tagRecord );

    // This puts characters from xmlStream into destinationWithoutHalt if it is
    // not NULL and destinationIncludingHalt if it is not NULL (if both are
    // NULL, the characters are discarded) until it reads in the first instance
    // of any of the characters in haltCharacters (which is put into
    // destinationIncludingHalt but not into destinationWithoutHalt).
    bool ReadToNextHalt( std::string const& haltCharacters,
                         std::ostream* destinationWithoutHalt,
                         std::ostream* destinationIncludingHalt );

    // This reads xmlStream until the end of the start or empty-element tag is
    // reached, parsing attributes along the way into attributeDestination,
    // throwing an exception if it does not reach the end of the tag without
    // the stream ending or finding malformed XML. All the characters of the
    // tag get put into tagRecord if it is not NULL.
    void CloseStartTag( AttributeMap* attributeDestination,
                        std::ostream* tagRecord );

    // This returns false if currentCharacter is '>' or if it is '/' and the
    // next character is '>', first setting currentElementIsEmpty
    // appropriately. It throws an exception if currentCharacter is '/' and is
    // not followed immediately by '>'. It returns true otherwise.
    bool TagIsStillOpen( std::ostream* tagRecord );

    // This reads characters from xmlStream into currentChar until the first
    // non-whitespace character is read in, returning true unless the stream
    // ended without any non-whitespace character being read in. If
    // destinationForReadCharacters is not NULL, all read characters (including
    // the first non-whitespace character) are put into it.
    bool SkipWhitespace( std::ostream* destinationForReadCharacters );

    // This parses an attribute from xmlStream assuming that the first
    // character of the attribute name is already in currentChar, and puts the
    // attribute into attributeDestination. It throws an exception if the XML
    // is malformed. All the characters of the tag read by this function get
    // put into tagRecord if it is not NULL.
    void ParseAttribute( AttributeMap* attributeDestination,
                         std::ostream* tagRecord );

    // This returns true if the attribute was correctly formed, assuming that
    // the first character of the attribute name is already in nameStream,
    // putting the rest of the name into nameStream and the value into
    // valueStream (without the quote marks), and all read characters into
    // tagRecord if it is not NULL.
    bool TryToReadValidAttribute( std::ostream& nameStream,
                                  std::ostream& valueStream,
                                  std::ostream* tagRecord )
    { return ( ReadToNextHalt( ( AllowedWhitespaceChars() + "=" ),
                               &nameStream,
                               tagRecord )
               &&
               SkipWhitespace( tagRecord )
               &&
               ( currentCharacter == '=' )
               &&
               ReadCharacter( tagRecord )
               &&
               SkipWhitespace( tagRecord )
               &&
               ( ( currentCharacter == '\'' ) || ( currentCharacter == '\"' ) )
               &&
               ReadToNextHalt( currentCharacter,
                               &valueStream,
                               tagRecord ) ); }

    // This puts everything except comments from xmlStream into currentBody up
    // to the end tag for the element named elementName. If any nested start
    // tags are found for child elements with this name, the text is read in
    // until each child of that name is closed and then to the next end tag for
    // that name. It does not validate that any other tags or elements are
    // correctly nested and so on. It throws an exception if it cannot find the
    // end of the element.
    void RecordToEndOfElement( std::string const& elementName );

    // This closes the markup just opened, assuming that currentCharacter is
    // the character just after '<'. It puts the content of the markup into
    // contentStream if it is not the end tag of the current element. The
    // markup is assumed to be the end tag of the current element if it is an
    // end tag for the element with name elementName and if
    // numberOfUnclosedElementsOfGivenName is 1. It returns the number
    // of elements with name given by elementName which have not yet had their
    // end tags found.
    unsigned int CloseMarkup( std::ostream& contentStream,
                              std::string const& elementName,
                      unsigned int const numberOfUnclosedElementsOfGivenName );

    // This reads in the markup just opened, assuming that currentCharacter is
    // the first character of a valid name of a start tag, putting the entire
    // tag into contentStream, and returning the number of elements with name
    // given by elementName which have not yet had their end tags found.
    unsigned int CloseStartTag( std::ostream& contentStream,
                                std::string const& elementName,
                      unsigned int const numberOfUnclosedElementsOfGivenName );

    // This reads in the markup just opened, assuming that currentCharacter is
    // '/' and that the markup is a valid end tag, putting the entire tag into
    // contentStream unless this was the end tag for the current open element,
    // and returning the number of elements with name given by elementName
    // which have not yet had their end tags found.
    unsigned int CloseEndTag( std::ostream& contentStream,
                                std::string const& elementName,
                      unsigned int const numberOfUnclosedElementsOfGivenName );

    // This tries to read an end tag, assuming that "</" was the last pair of
    // characters read in before this function was called, putting the name of
    // the element into nameDestination. The characters read in by this
    // function are passed to tagRecord. It throws an exception if it could not
    // read a valid end tag.
    void ReadEndTag( std::string& nameDestination,
                     std::ostream& tagRecord );
  };





  // This loads the given string as the content to parse into XML elements.
  inline void
  RestrictedXmlParser::LoadString( std::string const& stringToParse )
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
  inline void
  RestrictedXmlParser::OpenRootElementOfFile( std::string const& fileName )
  {
    CloseFile();
    xmlFileStream.open( fileName.c_str() );
    xmlStream = &xmlFileStream;
    ReadPrologAndOpenRootElement();
  }

  // This closes the file which was being parsed, if there is one open.
  inline void RestrictedXmlParser::CloseFile()
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
  inline void RestrictedXmlParser::ReadAllOfRootElementOfFile(
                                                  std::string const& fileName )
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
  inline bool RestrictedXmlParser::ReadNextElement()
  {
    if( ReadToNextTagOpener( NULL ) )
    {
      currentName.assign( "" );
      currentAttributes.clear();
      currentBody.assign( "" );
      currentElementIsEmpty = true;
      ReadStartTag( currentName,
                    &currentAttributes,
                    NULL );
      if( !currentElementIsEmpty )
      {
        RecordToEndOfElement( currentName );
      }
      return true;
    }
    return false;
  }

  // This resets the parser so that ReadNextElement() will start from the
  // beginning of the content to be parsed (which is either the string given
  // to LoadString(...) or the body of the root element of the file given to
  // OpenRootElementOfFile(...); if called after calling
  // ReadAllOfRootElementOfFile(...), nothing happens).
  inline void RestrictedXmlParser::ReturnToBeginningOfText()
  {
    ResetContent();
    xmlStream->clear();
    xmlStream->seekg( std::ios::beg );
  }

  // This sets the various recording data to the values they should have
  // before reading in some text.
  inline void RestrictedXmlParser::ResetContent()
  {
    fileProlog.assign( "" );
    rootName.assign( "" );
    rootAttributes.clear();
    currentName.assign( "" );
    currentAttributes.clear();
    currentBody.assign( "" );
    currentElementIsEmpty = true;
  }

  // This resets the content data members, then reads in everything up to
  // (but not including) the '<' of the first valid opening tag of an
  // element, then reads the name in that tag into rootName, and the
  // attributes into rootAttributes, and leaves xmlStream ready to read in
  // the first character after the '>' which closes the opening tag of the
  // root element.
  inline void RestrictedXmlParser::ReadPrologAndOpenRootElement()
  {
    ResetContent();
    std::stringstream prologStream;
    std::stringstream tagStream;
    if( !(ReadToNextTagOpener( &prologStream )) )
    {
      throw std::runtime_error( "No root element found in file!" );
    }
    fileProlog = prologStream.str();
    ReadStartTag( rootName,
                  &rootAttributes,
                  &tagStream );
  }

  // This puts characters from xmlStream into nonMarkupDestination if not
  // NULL, or discards the characters if nonMarkupDestination is NULL, until
  // it reads in the first instance of '<' followed by an allowed element
  // name start character. Neither the '<' nor the starting character are put
  // into nonMarkupDestination, but the starting character is left in
  // currentChar. True is returned if such a sequence is found, false if not
  // and no more characters can be read.
  inline bool RestrictedXmlParser::ReadToNextTagOpener(
                                           std::ostream* nonMarkupDestination )
  {
    while( ReadToNextHalt( '<',
                           nonMarkupDestination,
                           NULL )
           &&
           xmlStream->get( currentCharacter ).good() )
    {
      if( ParsingUtilities::CharacterIsInString( currentCharacter,
                                                 AllowedNameStartChars() ) )
      {
        return true;
      }
      else
      {
        TryToCloseNonTagMarkup( nonMarkupDestination );
      }
    }
    return false;
  }

  // This either reads in a processing instruction or a comment or a CDATA
  // section, passing relevant characters to contentStream, or just puts '<'
  // followed by currentCharacter into contentStream, depending on what
  // currentCharacter is.
  inline void RestrictedXmlParser::TryToCloseNonTagMarkup(
                                   std::ostream* destinationForReadCharacters )
  {
    if( currentCharacter == '?' )
    {
      CloseQuestionMark( destinationForReadCharacters );
    }
    else if( currentCharacter == '!' )
    {
      CloseExclamationMark( destinationForReadCharacters );
    }
    else
    {
      if( destinationForReadCharacters != NULL )
      {
        (*destinationForReadCharacters) << '<' << currentCharacter;
      }
    }
  }

  // This puts "<?" into destinationForReadCharacters if it is not NULL, then
  // reads in characters from xmlStream until the next "?>", placing
  // them all into destinationForReadCharacters if it is not NULL. It throws
  // an exception if it fails to read in "?>".
  inline void RestrictedXmlParser::CloseQuestionMark(
                                   std::ostream* destinationForReadCharacters )
  {
    if( destinationForReadCharacters != NULL )
    {
      (*destinationForReadCharacters) << "<?";
    }
    while( ReadToNextHalt( '?',
                           NULL,
                           destinationForReadCharacters )
           &&
           ReadCharacter( destinationForReadCharacters ) )
    {
      if( currentCharacter == '>' )
      {
        return;
      }
    }
    throw std::runtime_error( "Failed to close processing instruction!" );
  }

  // This puts characters from xmlStream into destinationWithoutHalt if it is
  // not NULL and destinationIncludingHalt if it is not NULL (if both are
  // NULL, the characters are discarded) until it reads in the first instance
  // of haltCharacter (which is put into destinationIncludingHalt but not
  // into destinationWithoutHalt).
  inline bool RestrictedXmlParser::ReadToNextHalt( char const haltCharacter,
                                          std::ostream* destinationWithoutHalt,
                                       std::ostream* destinationIncludingHalt )
  {
    while( ReadCharacter( destinationIncludingHalt )
           &&
           ( currentCharacter != haltCharacter ) )
    {
      if( destinationWithoutHalt != NULL )
      {
        (*destinationWithoutHalt) << currentCharacter;
      }
    }
    return ( currentCharacter == haltCharacter );
  }

  // This puts the next character from xmlStream into currentChar, then puts
  // it into destinationForReadCharacters if it is not NULL, then returns
  // true, unless no character could be read, in which case false is
  // returned.
  inline bool RestrictedXmlParser::ReadCharacter(
                                    std::ostream* destinationForReadCharacter )
  {
    if( xmlStream->get( currentCharacter ).good() )
    {
      if( destinationForReadCharacter != NULL )
      {
        (*destinationForReadCharacter) << currentCharacter;
      }
      return true;
    }
    else
    {
      return false;
    }
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
  inline void RestrictedXmlParser::CloseExclamationMark(
                                   std::ostream* destinationForReadCharacters )
  {
    if( !(xmlStream->get( currentCharacter ).good()) )
    {
      throw
      std::runtime_error( "Failed to close structure following \"<!\"" );
    }
    if( currentCharacter == '-' )
    {
      ReadComment();
    }
    else if( currentCharacter == '[' )
    {
      ReadCdata( destinationForReadCharacters );
    }
    else
    {
      if( destinationForReadCharacters != NULL )
      {
        (*destinationForReadCharacters) << '<' << currentCharacter;
      }
    }
  }

  // This tries to read in characters (discarding them) from xmlStream until
  // a valid comment has been read, assuming that the previous characters
  // were "<!-". It throws an exception if it fails to close the comment.
  inline void RestrictedXmlParser::ReadComment()
  {
    while( xmlStream->get( currentCharacter ).good() )
    {
      if( ( currentCharacter == '-' )
          &&
          xmlStream->get( currentCharacter ).good()
          &&
          ( currentCharacter == '-' )
          &&
          xmlStream->get( currentCharacter ).good()
          &&
          ( currentCharacter == '>' ) )
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
  inline void
  RestrictedXmlParser::ReadCdata( std::ostream* destinationForReadCharacters )
  {
    if( destinationForReadCharacters != NULL )
    {
      (*destinationForReadCharacters) << "<![";
    }
    while( ReadToNextHalt( ']',
                           NULL,
                           destinationForReadCharacters )
           &&
           ReadCharacter( destinationForReadCharacters )
           &&
           ( currentCharacter == ']' )
           &&
           ReadCharacter( destinationForReadCharacters ) )
    {
      if( currentCharacter == '>' )
      {
        return;
      }
    }
    throw std::runtime_error( "Failed to close <[...]]> structure!" );
  }

  // This tries to read a start tag, assuming that '<' was read before
  // currentChar, and that currentChar is the first character of a valid
  // element name. It sets currentElementIsEmpty to be true if the tag was
  // actually an empty-element tag rather than a start tag. All the
  // characters of the tag are passed into tagRecord if it is not NULL
  // (clearing the content of tagRecord beforehand). If it fails to read a
  // valid start tag, it throws an exception.
  inline void RestrictedXmlParser::ReadStartTag( std::string& nameDestination,
                                            AttributeMap* attributeDestination,
                                                 std::stringstream* tagRecord )
  {
    if( tagRecord != NULL )
    {
      ParsingUtilities::ResetStringstream( *tagRecord );
      (*tagRecord) << '<' << currentCharacter;
    }
    std::stringstream bufferStream;
    bufferStream << currentCharacter;

    // We read to the first character that marks the end of the element name.
    if( !(ReadToNextHalt( ( AllowedWhitespaceChars() + ">/" ),
                          &bufferStream,
                          tagRecord )) )
    {
      throw std::runtime_error(
             "Failed to find whitespace or end of tag after element name!" );
    }
    nameDestination = bufferStream.str();
    CloseStartTag( attributeDestination,
                   tagRecord );
  }

  // This puts characters from xmlStream into destinationWithoutHalt if it is
  // not NULL and destinationIncludingHalt if it is not NULL (if both are
  // NULL, the characters are discarded) until it reads in the first instance
  // of any of the characters in haltCharacters (which is put into
  // destinationIncludingHalt but not into destinationWithoutHalt).
  inline bool
  RestrictedXmlParser::ReadToNextHalt( std::string const& haltCharacters,
                                       std::ostream* destinationWithoutHalt,
                                       std::ostream* destinationIncludingHalt )
  {
    while( ReadCharacter( destinationIncludingHalt )
           &&
           !(ParsingUtilities::CharacterIsInString( currentCharacter,
                                                    haltCharacters )) )
    {
      if( destinationWithoutHalt != NULL )
      {
        (*destinationWithoutHalt) << currentCharacter;
      }
    }
    return ParsingUtilities::CharacterIsInString( currentCharacter,
                                                  haltCharacters );
  }

  // This reads xmlStream until the end of the start or empty-element tag is
  // reached, parsing attributes along the way into attributeDestination,
  // throwing an exception if it does not reach the end of the tag without
  // the stream ending or finding malformed XML. All the characters of the
  // tag get put into tagRecord if it is not NULL.
  inline void
  RestrictedXmlParser::CloseStartTag( AttributeMap* attributeDestination,
                                      std::ostream* tagRecord )
  {
    if( !(SkipWhitespace( tagRecord )) )
    {
      throw std::runtime_error( "Could not find end of start tag!" );
    }
    // Now currentChar is the end of the start tag ('>'), or the first
    // character of the end of an empty element tag ('/' in "/>") (or the
    // tag is malformed with a '/' out of place), or is the first character
    // of an attribute name.
    if( TagIsStillOpen( tagRecord ) )
    {
      ParseAttribute( attributeDestination,
                      tagRecord );
      if( !(ReadCharacter( tagRecord )) )
      {
        throw std::runtime_error( "Could not find end of start tag!" );
      }
      CloseStartTag( attributeDestination,
                     tagRecord );
    }
  }

  // This returns false if currentCharacter is '>' or if it is '/' and the
  // next character is '>', first setting currentElementIsEmpty
  // appropriately. It throws an exception if currentCharacter is '/' and is
  // not followed immediately by '>'. It returns true otherwise.
  inline bool
  RestrictedXmlParser::TagIsStillOpen( std::ostream* tagRecord )
  {
    if( currentCharacter == '>' )
    {
      currentElementIsEmpty = false;
      return false;
    }
    else if( currentCharacter == '/' )
    {
      if( ReadCharacter( tagRecord )
          &&
          ( currentCharacter == '>' ) )
      {
        currentElementIsEmpty = true;
        return false;
      }
      else
      {
        throw std::runtime_error( "Attribute name cannot begin with \'/\'" );
      }
    }
    return true;
  }

  // This reads characters from xmlStream into currentChar until the first
  // non-whitespace character is read in, returning true unless the stream
  // ended without any non-whitespace character being read in. If
  // destinationForReadCharacters is not NULL, all read characters (including
  // the first non-whitespace character) are put into it.
  inline bool RestrictedXmlParser::SkipWhitespace(
                                   std::ostream* destinationForReadCharacters )
  {
    // If the name is followed by whitespace, we keep going to the first
    // non-whitespace character, discarding all the whitespace characters
    // along the way.
    while( ParsingUtilities::CharacterIsInString( currentCharacter,
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
  // is malformed. All the characters of the tag read by this function get
  // put into tagRecord if it is not NULL.
  inline void
  RestrictedXmlParser::ParseAttribute( AttributeMap* attributeDestination,
                                       std::ostream* tagRecord )
  {
    std::stringstream nameStream;
    nameStream << currentCharacter;
    std::stringstream valueStream;
    if( !( TryToReadValidAttribute( nameStream,
                                    valueStream,
                                    tagRecord ) ) )
    {
      throw std::runtime_error( "Could not parse an attribute correctly!" );
    }
    if( attributeDestination != NULL )
    {
      (*attributeDestination)[ nameStream.str() ] = valueStream.str();
    }
  }


  // This puts everything except comments from xmlStream into currentBody up
  // to the end tag for the element named elementName. If any nested start
  // tags are found for child elements with this name, the text is read in
  // until each child of that name is closed and then to the next end tag for
  // that name. It does not validate that any other tags or elements are
  // correctly nested and so on. It throws an exception if it cannot find the
  // end of the element.
  inline void
  RestrictedXmlParser::RecordToEndOfElement( std::string const& elementName )
  {
    unsigned int numberOfUnclosedElementsOfGivenName( 1 );
    std::stringstream contentStream;
    while( ( numberOfUnclosedElementsOfGivenName > 0 )
           &&
           ReadToNextHalt( '<',
                           &contentStream,
                           NULL ) )
    {
      if( !(xmlStream->get( currentCharacter ).good()) )
      {
        std::stringstream errorBuilder;
        errorBuilder
        << "Could not find end tag for element <" << elementName << ">";
        throw std::runtime_error( errorBuilder.str() );
      }
      numberOfUnclosedElementsOfGivenName = CloseMarkup( contentStream,
                                                         elementName,
                                       numberOfUnclosedElementsOfGivenName );
    }
    currentBody = contentStream.str();
  }

  // This closes the markup just opened, assuming that currentCharacter is
  // the character just after '<'. It puts the content of the markup into
  // contentStream if it is not the end tag of the current element. The
  // markup is assumed to be the end tag of the current element if it is an
  // end tag for the element with name elementName and if
  // numberOfUnclosedElementsOfGivenName is 1. It returns the number
  // of elements with name given by elementName which have not yet had their
  // end tags found.
  inline unsigned int
  RestrictedXmlParser::CloseMarkup( std::ostream& contentStream,
                                    std::string const& elementName,
                       unsigned int const numberOfUnclosedElementsOfGivenName )
  {
    if( ParsingUtilities::CharacterIsInString( currentCharacter,
                                                  AllowedNameStartChars() ) )
    {
      return CloseStartTag( contentStream,
                            elementName,
                            numberOfUnclosedElementsOfGivenName );
    }
    else if( currentCharacter == '/' )
    {
      return CloseEndTag( contentStream,
                          elementName,
                          numberOfUnclosedElementsOfGivenName );
    }
    else
    {
      TryToCloseNonTagMarkup( &contentStream );
      return numberOfUnclosedElementsOfGivenName;
    }
  }

  // This reads in the markup just opened, assuming that currentCharacter is
  // the first character of a valid name of a start tag, putting the entire
  // tag into contentStream, and returning the number of elements with name
  // given by elementName which have not yet had their end tags found.
  inline unsigned int
  RestrictedXmlParser::CloseStartTag( std::ostream& contentStream,
                                      std::string const& elementName,
                       unsigned int const numberOfUnclosedElementsOfGivenName )
  {
    std::stringstream tagStream;
    std::string innerName;
    ReadStartTag( innerName,
                  NULL,
                  &tagStream );
    contentStream << tagStream.str();
    return ( ( !currentElementIsEmpty && ( innerName == elementName ) ) ?
             ( numberOfUnclosedElementsOfGivenName + 1 ):
             numberOfUnclosedElementsOfGivenName );
  }

  // This reads in the markup just opened, assuming that currentCharacter is
  // '/' and that the markup is a valid end tag, putting the entire tag into
  // contentStream unless this was the end tag for the current open element,
  // and returning the number of elements with name given by elementName
  // which have not yet had their end tags found.
  inline unsigned int
  RestrictedXmlParser::CloseEndTag( std::ostream& contentStream,
                                    std::string const& elementName,
                       unsigned int const numberOfUnclosedElementsOfGivenName )
  {
    std::stringstream tagStream;
    std::string innerName;
    ReadEndTag( innerName,
                tagStream );
    if( innerName == elementName )
    {
      if( numberOfUnclosedElementsOfGivenName > 1 )
      {
        contentStream << tagStream.str();
      }
      return ( numberOfUnclosedElementsOfGivenName - 1 );
    }
    else
    {
      contentStream << tagStream.str();
      return numberOfUnclosedElementsOfGivenName;
    }
  }

  // This tries to read an end tag, assuming that "</" was the last pair of
  // characters read in before this function was called, putting the name of
  // the element into nameDestination. The characters read in by this
  // function are passed to tagRecord. It throws an exception if it could not
  // read a valid end tag.
  inline void RestrictedXmlParser::ReadEndTag( std::string& nameDestination,
                                               std::ostream& tagRecord )
  {
    tagRecord << "</";
    std::stringstream nameStream;
    if( !( ReadToNextHalt( ( AllowedWhitespaceChars() + ">" ),
                           &nameStream,
                           &tagRecord )
           &&
           SkipWhitespace( &tagRecord )
           &&
           ( currentCharacter == '>' ) ) )
    {
      throw std::runtime_error( "Could not close a valid end tag!" );
    }
    nameDestination = nameStream.str();
  }

}

#endif /* RESTRICTEDXMLPARSER_HPP_ */
