/*
 * AsciiXmlParser.hpp
 *
 *  Created on: Dec 8, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef ASCIIXMLPARSER_HPP_
#define ASCIIXMLPARSER_HPP_

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
  //
  // This class is compatible with UTF-8 in a sense, as the bytes will be read
  // correctly, but they will be stored in string objects, which are
  // effectively arrays of char data, so (depending on architecture in a way
  // that is beyond my expertise) will be stored as a series of bytes.
  // Conversion back to UTF-8 is beyond the scope of this class. The Boost
  // library provides string UTF conversion tools.
  // (http://www.boost.org/doc/libs/1_49_0/libs/locale/doc/html/charset_handling.html)
  //
  // This references http://www.w3.org/TR/2008/REC-xml-20081126/ for XML
  // definition. It does not comply with the requirement that it MUST accept
  // UTF-8 and UTF-16, and does not do a huge amount of validation.
  //
  // It's worth checking this out:
  // http://www.joelonsoftware.com/articles/Unicode.html
  //
  // Notes to self:
  // #x9 is '\t'
  // #xA is '\n'
  // #xD is '\r'
  // #x20 is ' '
  class AsciiXmlParser
  {
  public:
    AsciiXmlParser( bool const isVerbose = false );
    ~AsciiXmlParser();


    // This loads the given string as the content to parse into XML elements.
    void LoadString( std::string const& stringToParse )
    {
      ResetContent();
      xmlStringStream.clear();
      xmlStringStream.str( stringToParse );
      xmlStream = &xmlStringStream;
    }

    // This opens the file with name fileName, and reads until the first
    // element opening tag is found, storing the prolog, and reads in that
    // first opening tag, then stops. The prolog is accessible through
    // FileProlog(), and the content of the opening tag is accessible through
    // RootName() and RootAttributes(). The rest of the root element is treated
    // as the content to parse into XML elements, as stringToParse would be in
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
      RecordToEndOfElement();
    }

    // This moves forward in the current content to parse until the next
    // element opening tag is found, then that element is read in until its
    // closing tag is found, and the character just after that is where the
    // next call of this function will resume. The element just parsed is
    // accessible through CurrentName() for its name, CurrentAttributes() for
    // its attributes, and CurrentBody() for the text between the opening and
    // closing tags. The return value is true if an element was read in, false
    // if no more valid elements could be found. (If no further elements are
    // found, the next call of this function will do nothing, as the next
    // character will just be the EOF character.)
    bool ReadNextElement() { return ( DiscardToNextTag()
                                      &&
                                      ParseTagName( currentName )
                                      &&
                                      ParseAttributes( currentAttributes )
                                      &&
                                      RecordToEndOfElement() ); }

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

    // This returns the name of the root element of the file opened by
    // OpenRootElementOfFile(...) or ReadAllOfRootElementOfFile(...).
    std::string const& CurrentName() { return currentName; }

    // This returns the attributes of the root element of the file opened by
    // OpenRootElementOfFile(...) or ReadAllOfRootElementOfFile(...).
    std::map< std::string, std::string > const& CurrentAttributes()
    { return currentAttributes; }

    // This returns the name of the root element of the file opened by
    // OpenRootElementOfFile(...) or ReadAllOfRootElementOfFile(...).
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
    static char MarkupOpener() { return '<'; }
    static char const MarkupCloser() { return '>'; }
    static char const TagCloser() { return '/'; }
    static std::string const AllowedWhitespaceChars() { return " \t\r\n"; }
    static std::string const AllowedQuoteChars() { return "\'\""; }
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
    std::stringstream bufferStream;
    char currentChar;

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
      while( ReadToNextTagOpener( &prologStream )
             &&
             !(TryToReadOpeningTag( rootName,
                                    &rootAttributes )) )
      {
        prologStream << bufferStream.str();
      }
    }

    // This puts characters from xmlStream into nonMarkupDestination if not
    // NULL, or discards the characters if nonMarkupDestination is NULL, until
    // it reads in the first instance of haltCharacter (which is not put into
    // nonMarkupDestination even if nonMarkupDestination is not NULL).
    bool ReadToNextTagOpener( std::istream* nonMarkupDestination )
    {
      while( ReadToNextHalt( '<',
                             nonMarkupDestination ) )
      {
        if( xmlStream->get( currentChar ).good() )
        {
          if( currentChar == '?' )
          {
            CloseQuestionMark( nonMarkupDestination );
          }
          else if( currentChar == '!' )
          {
            CloseExclamationMark( nonMarkupDestination );
          }
          else if( !(ParsingUtilities::CharacterIsInString( currentChar,
                                                     ( AllowedWhitespaceChars()
                                                   + AllowedQuoteChars() ) )) )
          {
            return true;
          }

        }
        else
        {
          return false;
        }
      }
      return false;
    }

    //
    void CloseQuestionMark( std::istream* nonMarkupDestination )
    {
      if( nonMarkupDestination != NULL )
      {
        nonMarkupDestination << "<?";
      }
      while( ReadToNextHalt( '?',
                             nonMarkupDestination ) )
      {
        if( xmlStream->get( currentChar ).good() )
        {
          if( nonMarkupDestination != NULL )
          {
            nonMarkupDestination << '?' << currentChar;
          }
          if( currentChar == '>' )
          {
            return;
          }
        }
      }
    }

    //
    void CloseExclamationMark( std::istream* nonMarkupDestination )
    {
      bufferStream.clear();
      bufferStream.str( "<!" );
      if( xmlStream->get( currentChar ).good() )
      {
        bufferStream << currentChar;
        if( ( currentChar == '-' )
            &&
            xmlStream->get( currentChar ).good() )
        {
          if( currentChar == '-' )
          {
            CloseComment();
          }
          else
          {
            nonMarkupDestination << "<!-" << currentChar;
          }
        }
        else if( ( currentChar == '[' )
                 &&
                 xmlStream->get( currentChar ).good() )
        {
          if( ( currentChar == 'C' )
              &&
              xmlStream->get( currentChar ).good() )
          {
            if( ( currentChar == 'D' )
                &&
                xmlStream->get( currentChar ).good() )
            {
              if( ( currentChar == 'A' )
                  &&
                  xmlStream->get( currentChar ).good() )
              {
                if( ( currentChar == 'T' )
                    &&
                    xmlStream->get( currentChar ).good() )
                {
                  if( ( currentChar == 'A' )
                      &&
                      xmlStream->get( currentChar ).good() )
                  {
                    if( ( currentChar == '[' )
                        &&
                        xmlStream->get( currentChar ).good() )
                    {
                      CloseCdata();
                    }
                    else
                    {
                      nonMarkupDestination << "<![CDATA" << currentChar;
                    }
                  }
                  else
                  {
                    nonMarkupDestination << "<![CDAT" << currentChar;
                  }
                }
                else
                {
                  nonMarkupDestination << "<![CDAT" << currentChar;
                }
              }
              else
              {
                nonMarkupDestination << "<![CDA" << currentChar;
              }
            }
            else
            {
              nonMarkupDestination << "<![CD" << currentChar;
            }
          }
          else
          {
            nonMarkupDestination << "<![C" << currentChar;
          }
        }
        else
        {
          nonMarkupDestination << "<![" << currentChar;
        }
      }
    }

    // This puts characters from xmlStream into nonMarkupDestination if not
    // NULL, or discards the characters if nonMarkupDestination is NULL, until
    // it reads in the first instance of haltCharacter (which is not put into
    // nonMarkupDestination even if nonMarkupDestination is not NULL).
    bool ReadToNextHalt( char const haltCharacter,
                         std::istream* nonMarkupDestination )
    {
      while( xmlStream->get( currentChar ).good()
             &&
             currentChar != haltCharacter )
      {
        if( nonMarkupDestination != NULL )
        {
          (*nonMarkupDestination) << currentChar;
        }
      }
      return ( currentChar == haltCharacter );
    }

    //
    bool TryToReadOpeningTag( std::string& nameDestination,
                              AttributeMap* attributeDestination )
    {
      // currentChar is '<' and xmlStream->get( currentChar ) will put the next
      // character into it.
      xmlStream->get( currentChar );
      bufferStream.clear();
      bufferStream.str( std::string( 1,
                                     currentChar ) );
      // There are a lot of characters which are specified as allowed initial
      // characters for an element name, but it's easier here to allow plenty
      // of technically invalid names, while deciding what is not an element
      // based on the first character after the '<'.
      if( ParsingUtilities::CharacterIsInString( currentChar,
                                                 "? \n\r\t" ) )
      {
        return false;
      }
      else if( currentChar == '!' )
      {
        DiscardPossibleComment();
        return false;
      }

      // At this point, currentChar is not whitespace or a character which
      // indicates the start of a comment or CDATA or process instruction etc.,
      // so it is taken as the start of the name.
    }

    // This puts characters from xmlStream into nonMarkupDestination if not
    // NULL, or discards the characters if nonMarkupDestination is NULL, until
    // it reads in the first instance of any of the characters in
    // haltCharacters (which is not put into nonMarkupDestination even if
    // nonMarkupDestination is not NULL).
    bool ReadToNext( std::string const& haltCharacters,
                     std::istream* nonMarkupDestination )
    {
      while( xmlStream->get( currentChar ).good()
             &&
             !(ParsingUtilities::CharacterIsInString( currentChar,
                                                      haltCharacters )) )
      {
        if( nonMarkupDestination != NULL )
        {
          (*nonMarkupDestination) << currentChar;
        }
      }
      return ParsingUtilities::CharacterIsInString( currentChar,
                                                    haltCharacters );
    }

    //
    void RecordToEndOfElement();

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

#endif /* ASCIIXMLPARSER_HPP_ */
