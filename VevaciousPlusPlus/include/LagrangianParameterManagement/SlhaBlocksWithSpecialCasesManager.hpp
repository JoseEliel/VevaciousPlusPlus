/*
 * SlhaBlocksWithSpecialCasesManager.hpp
 *
 *  Created on: Nov 2, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHABLOCKSWITHSPECIALCASESMANAGER_HPP_
#define SLHABLOCKSWITHSPECIALCASESMANAGER_HPP_

#include "CommonIncludes.hpp"
#include "LesHouchesAccordBlockEntryManager.hpp"
#include "SlhaSourcedParameterFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaDsbHiggsVevFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaHiggsMixingBilinearFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaMassSquaredDiagonalFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaTrilinearDiagonalFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaBlocksWithSpecialCasesManager :
                                       public LesHouchesAccordBlockEntryManager
  {
  public:
    SlhaBlocksWithSpecialCasesManager( std::string const& validBlocksString,
                                       std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                       std::string const& fixedScaleType,
                                       std::string const& fixedScaleArgument );
    SlhaBlocksWithSpecialCasesManager(
                                 std::set< std::string > const& validBlocksSet,
                                       std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                       std::string const& fixedScaleType,
                                       std::string const& fixedScaleArgument );
    virtual ~SlhaBlocksWithSpecialCasesManager();


  protected:
    std::vector< SlhaSourcedParameterFunctionoid* > activeDerivedParameters;
    std::map< std::string, std::string > aliasesToSwitchStrings;

    // This adds all the valid aliases to aliasesToSwitchStrings.
    void InitializeSlhaOneOrTwoAliases();

    // This adds newParameter to activeDerivedParameters and updates
    // numberOfDistinctActiveParameters and activeParametersToIndices (all
    // aliases in aliasesToSwitchStrings which map to switchString are added
    // mapping to the index of this parameter), and returns true paired with
    // the index of newParameter. (The index mapped to in
    // activeParametersToIndices and used as part of the return value is based
    // on numberOfDistinctActiveParameters before it gets incremented, rather
    // than using newParameter->IndexInValuesVector().)
    std::pair< bool, size_t >
    AddNewDerivedParameter( std::string const& switchString,
                            SlhaSourcedParameterFunctionoid* newParameter );

    // This checks parameterName against all the special cases. If the
    // parameter name is not recognized as a valid special case, then the
    // returned pair is false paired with (size_t)(-1).
    virtual std::pair< bool, size_t >
    RegisterUnregisteredParameter( std::string const& parameterName );

    // This adds the parameter based on the alias given by switchString for the
    // parameter.
    std::pair< bool, size_t >
    SlhaOneOrTwoSpecialCase( std::string const& switchString );

    // This registers the required TX[0,0], AX[0,0], and YX[0,0] parameters
    // where X is replaced by sfermionType and 0 by generationIndex, and then
    // returns the result of AddNewDerivedParameter using a new
    // SlhaTrilinearDiagonalFunctionoid constructed with those parameter
    // functionoids.
    std::pair< bool, size_t >
    RegisterSlhaOneOrTwoCompatibleTrilinear( std::string const& switchString,
                                             char const sfermionType,
                                             char const generationIndex );

    // This registers the required M2X[0,0] and MSOFT[msoftIndex] parameters
    // where X is replaced by sfermionType and 0 by generationIndex, and then
    // returns the result of AddNewDerivedParameter using a new
    // SlhaMassSquaredDiagonalFunctionoid constructed with those parameter
    // functionoids.
    std::pair< bool, size_t >
    RegisterSlhaOneOrTwoCompatibleMassSquared( std::string const& switchString,
                                               char const sfermionType,
                                               char const generationIndex,
                                               std::string const& msoftIndex );

    // This maps a set of strings to switchString in aliasesToSwitchStrings:
    // switchString itself, and every upper-/lowercase combination of
    // blockString (after the indices bracket has been put in the correct
    // format, and only changing the cases of the block name characters).
    void CoverAllCasesForSlhaBlock( std::string const& switchString,
                                    std::string const& blockString );

    // This finds the functionoid in activeDerivedParameters which has the
    // given index (which is not its position in the activeDerivedParameters
    // vector, importantly).
    SlhaSourcedParameterFunctionoid const*
    FindActiveDerivedParameter( size_t const soughtIndex ) const;
  };

  // This adds newParameter to activeDerivedParameters and updates
  // numberOfDistinctActiveParameters and activeParametersToIndices (all
  // aliases in aliasesToSwitchStrings which map to switchString are added
  // mapping to the index of this parameter), and returns true paired with
  // the index of newParameter. (The index mapped to in
  // activeParametersToIndices and used as part of the return value is based
  // on numberOfDistinctActiveParameters before it gets incremented, rather
  // than using newParameter->IndexInValuesVector().)
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::AddNewDerivedParameter(
                                               std::string const& switchString,
                                SlhaSourcedParameterFunctionoid* newParameter )
  {
    activeDerivedParameters.push_back( newParameter );
    for( std::map< std::string, std::string >::const_iterator
         aliasToSwitchString( aliasesToSwitchStrings.begin() );
         aliasToSwitchString < aliasesToSwitchStrings.end();
         ++aliasToSwitchString )
    {
      if( aliasToSwitchString->second == switchString )
      {
        activeParametersToIndices[ aliasToSwitchString->first ]
        = numberOfDistinctActiveParameters;
      }
    }
    // We tersely update numberOfDistinctActiveParameters with
    // post-increment ++ so that the pair has the correct index.
    return std::pair< bool, size_t >( true,
                                      numberOfDistinctActiveParameters++ );
  }

  // This checks parameterName against all the special cases. If the
  // parameter name is not recognized as a valid special case, then the
  // returned pair is false paired with (size_t)(-1).
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::RegisterUnregisteredParameter(
                                             std::string const& parameterName )
  {
    std::map< std::string, std::string >::const_iterator
    aliasToSwitchString( aliasesToSwitchStrings.find( parameterName ) );
    if( aliasToSwitchString == aliasesToSwitchStrings.end() )
    {
      return std::pair< bool, size_t >( false,
                                        -1 );
    }
    return SlhaOneOrTwoSpecialCase( aliasToSwitchString->second );
  }


  // This registers the required TX[0,0], AX[0,0], and YX[0,0] parameters
  // where X is replaced by sfermionType and 0 by generationIndex, and then
  // returns the result of AddNewDerivedParameter using a new
  // SlhaTrilinearDiagonalFunctionoid constructed with those parameter
  // functionoids.
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::RegisterSlhaOneOrTwoCompatibleTrilinear(
                                               std::string const& switchString,
                                                       char const sfermionType,
                                                   char const generationIndex )
  {
    std::string blockNameWithIndices( "TX[0,0]" );
    blockNameWithIndices[ 1 ] = sfermionType;
    blockNameWithIndices[ 3 ] = blockNameWithIndices[ 5 ] = generationIndex;
    SlhaSourcedParameterFunctionoid const&
    directTrilinear( RegisterBlockEntry( blockNameWithIndices ) );
    blockNameWithIndices[ 0 ] = 'A';
    SlhaSourcedParameterFunctionoid const&
    trilinearOverYukawa( RegisterBlockEntry( blockNameWithIndices ) );
    blockNameWithIndices[ 0 ] = 'Y';
    SlhaSourcedParameterFunctionoid const&
    appropriateYukawa( RegisterBlockEntry( blockNameWithIndices ) );
    return AddNewDerivedParameter( switchString,
                                   new SlhaTrilinearDiagonalFunctionoid(
                                            numberOfDistinctActiveParameters,
                                                             directTrilinear,
                                                         trilinearOverYukawa,
                                                       appropriateYukawa ) );
  }

  // This registers the required M2X[0,0] and MSOFT[msoftIndex] parameters
  // where X is replaced by sfermionType and 0 by generationIndex, and then
  // returns the result of AddNewDerivedParameter using a new
  // SlhaMassSquaredDiagonalFunctionoid constructed with those parameter
  // functionoids.
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::RegisterSlhaOneOrTwoCompatibleMassSquared(
                                               std::string const& switchString,
                                                       char const sfermionType,
                                                    char const generationIndex,
                                                std::string const& msoftIndex )
  {
    std::string blockNameWithIndices( "MSX2[0,0]" );
    blockNameWithIndices[ 2 ] = sfermionType;
    blockNameWithIndices[ 5 ] = blockNameWithIndices[ 7 ] = generationIndex;
    SlhaSourcedParameterFunctionoid const&
    squareMass( RegisterBlockEntry( blockNameWithIndices ) );
    blockNameWithIndices.assign( "MSOFT[" );
    blockNameWithIndices.append( msoftIndex );
    blockNameWithIndices.append( "]" );
    SlhaSourcedParameterFunctionoid const&
    linearMass( RegisterBlockEntry( blockNameWithIndices  ) );
    return AddNewDerivedParameter( switchString,
                                   new SlhaMassSquaredDiagonalFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                                    squareMass,
                                                                linearMass ) );
  }

  // This maps a set of strings to switchString in aliasesToSwitchStrings:
  // switchString itself, and every upper-/lowercase combination of
  // blockString (after the indices bracket has been put in the correct
  // format, and only changing the cases of the block name characters).
  inline void SlhaBlocksWithSpecialCasesManager::CoverAllCasesForSlhaBlock(
                                               std::string const& switchString,
                                               std::string const& blockString )
  {
    aliasesToSwitchStrings[ switchString ] = switchString;
    std::string blockAfterFormat( FormatVariable( blockString ) );
    size_t const lastBlockNameCharacter( blockAfterFormat.find( '[',
                                                                1 ) );
    for( size_t characterIndex( 0 );
         characterIndex < lastBlockNameCharacter;
         ++characterIndex )
    {
      blockAfterFormat[ characterIndex ]
      = toupper( blockAfterFormat[ characterIndex ] );
      aliasesToSwitchStrings[ blockAfterFormat ] = switchString;
      blockAfterFormat[ characterIndex ]
      = tolower( blockAfterFormat[ characterIndex ] );
      aliasesToSwitchStrings[ blockAfterFormat ] = switchString;
    }
  }

  // This finds the functionoid in activeDerivedParameters which has the
  // given index (which is not its position in the activeDerivedParameters
  // vector, importantly).
  inline SlhaSourcedParameterFunctionoid const*
  SlhaBlocksWithSpecialCasesManager::FindActiveDerivedParameter(
                                               size_t const soughtIndex ) const
  {
    // We unfortunately need to do a linear search through
    // activeDerivedParameters to find the functionoid with the matching
    // index.
    for( std::vector< SlhaSourcedParameterFunctionoid* >::const_iterator
         activeParameter( activeDerivedParameters.begin() );
         activeParameter < activeDerivedParameters.end();
         ++activeParameter )
    {
      if( (*activeParameter)->IndexInValuesVector() == soughtIndex )
      {
        return (*activeParameter);
      }
    }
    std::stringstream errorBuilder;
    errorBuilder
    << "SlhaBlocksWithSpecialCasesManager::FindActiveDerivedParameter( "
    << soughtIndex << " ) found no SlhaSourcedParameterFunctionoid* in"
    << " activeDerivedParameters with that index.";
    throw std::out_of_range( errorBuilder.str() );
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHABLOCKSWITHSPECIALCASESMANAGER_HPP_ */
