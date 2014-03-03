/*
 * RunningParameterManager.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef RUNNINGPARAMETERMANAGER_HPP_
#define RUNNINGPARAMETERMANAGER_HPP_

#include "../StandardIncludes.hpp"
#include "ParameterFunctionoid.hpp"
#include "SLHA.hpp"
#include "SlhaFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class RunningParameterManager
  {
  public:
    typedef LHPC::SLHA::SparseManyIndexedBlock< double > RestrictedSlhaBlock;

    // This puts all variables with index brackets into a consistent form.
    static std::string
    FormatVariable( std::string const& unformattedVariable );


    RunningParameterManager();
    virtual
    ~RunningParameterManager();


    // This adds to the list of SLHA blocks that will be read in when updating
    // running parameters. An alias for the block can be provided.
    void AddValidSlhaBlock( std::string slhaBlock,
                            std::string const& aliasString );

    // This creates a new ParameterFunctionoid and a string mapping to it.
    void CreateDerivedParameter( std::string const& parameterString );

    // This returns a pointer to the ParameterFunctionoid mapped to by
    // parameterString, returning NULL if there is none.
    ParameterFunctionoid* GetFunctionoid( std::string const& parameterString );

    // This updates all SlhaFunctionoids with data from a new SLHA file.
    void UpdateSlhaParameters( std::string const& slhaFilename );

    // This updates all functionoids in order based on the new renormalization
    // scale.
    void UpdateRunningParameters( double const renormalizationScale );


  protected:

    // This puts all index brackets into a consistent form.
    static std::string
    FormatIndexBracketContent( std::string const& unformattedBracketContent );

    std::vector< ParameterFunctionoid* > parameterFunctionoidPointers;
    std::vector< SlhaFunctionoid* > slhaFunctionoidPointers;
    std::map< std::string, ParameterFunctionoid* > parameterFunctionoidMap;
    std::map< std::string, ParameterFunctionoid* >::iterator functionoidFinder;
    LHPC::SlhaParser slhaParser;
    std::vector< RestrictedSlhaBlock* > slhaBlockPointers;
    std::map< std::string, RestrictedSlhaBlock* > slhaBlockMap;
    std::map< std::string, RestrictedSlhaBlock* >::iterator slhaBlockFinder;
    std::map< std::string, std::string > slhaAliasMap;
    std::map< std::string, std::string >::iterator slhaAliasFinder;

    // This returns a pointer to the ParameterFunctionoid for the block named
    // blockName with indices given by indexString, making a new one if there
    // is none already. If the block name was invalid though, nothing is made
    // and NULL is returned.
    ParameterFunctionoid* GetSlhaFunctionoid( std::string const& blockName,
                                              std::string const& indexString );
  };




  // This puts all variables with index brackets into a consistent form.
  inline std::string RunningParameterManager::FormatVariable(
                                       std::string const& unformattedVariable )
  {
    size_t openBracket( unformattedVariable.find( '[' ) );
    if( openBracket == std::string::npos )
    {
      return std::string( unformattedVariable );
    }
    if( unformattedVariable[ unformattedVariable.size() - 1 ] != ']' )
    {
      throw std::runtime_error(
                         "In parsing model file, [...] not closed properly." );
    }
    std::stringstream formattedStream;
    formattedStream << unformattedVariable.substr( 0,
                                                   openBracket ) << '['
    << FormatIndexBracketContent( unformattedVariable.substr(
                                                           ( openBracket + 1 ),
                   ( unformattedVariable.size() - openBracket - 2 ) ) ) << ']';
    return std::string( formattedStream.str() );
  }

  // This returns a pointer to the ParameterFunctionoid mapped to by
  // parameterString, returning NULL if there is none.
  inline ParameterFunctionoid*
  RunningParameterManager::GetFunctionoid( std::string const& parameterString )
  {
    functionoidFinder = parameterFunctionoidMap.find( parameterString );
    if( functionoidFinder != parameterFunctionoidMap.end() )
    {
      return functionoidFinder->second;
    }
    size_t bracketPosition( parameterString.find_first_of( '[' ) );
    if( bracketPosition < ( parameterString.size() - 2 ) )
    {
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "parameterString = \"" << parameterString << "\"";
      std::cout << std::endl;/**/

      return GetSlhaFunctionoid( parameterString.substr( 0,
                                                         bracketPosition ),
                               parameterString.substr( ( bracketPosition + 1 ),
                          ( parameterString.size() - bracketPosition - 2 ) ) );
    }
    return NULL;
  }

  // This updates all SlhaFunctionoids with data from a new SLHA file.
  inline void RunningParameterManager::UpdateSlhaParameters(
                                              std::string const& slhaFilename )
  {
    slhaParser.readFile( slhaFilename );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "slhaBlockPointers.size() = " << slhaBlockPointers.size();
    for( std::vector< RestrictedSlhaBlock* >::iterator
         whichBlock( slhaBlockPointers.begin() );
         whichBlock < slhaBlockPointers.end();
         ++whichBlock )
    {
      std::cout
      << std::endl << "Block " << (*whichBlock)->getName() << std::endl;
      for( int whichScale( 1 );
           whichScale <= (*whichBlock)->getNumberOfCopiesWithDifferentScale();
           ++whichScale )
      {
        std::cout
        << "scale copy " << whichScale << " has scale "
        << (*(*whichBlock))[ whichScale ].getScale() << std::endl
        << "as string =" << std::endl
        << (*(*whichBlock))[ whichScale ].getAsString();
      }
    }
    std::cout << std::endl;/**/

    for( std::vector< SlhaFunctionoid* >::iterator
         whichParameter( slhaFunctionoidPointers.begin() );
         whichParameter < slhaFunctionoidPointers.end();
         ++whichParameter )
    {
      (*whichParameter)->UpdateForNewSlhaParameters();
    }
  }

  // This updates all functionoids in order based on the new renormalization
  // scale.
  inline void RunningParameterManager::UpdateRunningParameters(
                                            double const renormalizationScale )
  {
    double const logarithmOfScale( log( renormalizationScale ) );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "logarithmOfScale = " << logarithmOfScale;
    std::cout << std::endl;/**/

    for( std::vector< ParameterFunctionoid* >::iterator
         whichParameter( parameterFunctionoidPointers.begin() );
         whichParameter < parameterFunctionoidPointers.end();
         ++whichParameter )
    {
      (*whichParameter)->UpdateForNewLogarithmOfScale( logarithmOfScale );
    }
  }


  // This puts all index brackets into a consistent form.
  inline std::string RunningParameterManager::FormatIndexBracketContent(
                                 std::string const& unformattedBracketContent )
  {
    std::vector< int > indicesVector( BOL::StringParser::stringToIntVector(
                                                 unformattedBracketContent ) );
    std::stringstream indicesStream;
    for( std::vector< int >::iterator
         whichIndex( indicesVector.begin() );
         whichIndex < indicesVector.end();
         ++whichIndex )
    {
      if( whichIndex != indicesVector.begin() )
      {
        indicesStream << ',';
      }
      indicesStream << *whichIndex;
    }
    return std::string( indicesStream.str() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* RUNNINGPARAMETERMANAGER_HPP_ */
