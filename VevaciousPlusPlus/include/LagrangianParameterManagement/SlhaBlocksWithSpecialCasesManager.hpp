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


    // This is mainly for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    static bool
    SortParameterByIndex( SlhaSourcedParameterFunctionoid const* firstPointer,
                         SlhaSourcedParameterFunctionoid const* secondPointer )
    { return ( firstPointer->IndexInValuesVector()
               < secondPointer->IndexInValuesVector() ); }

    std::vector< SlhaSourcedParameterFunctionoid* > activeDerivedParameters;
    std::map< std::string, std::string > aliasesToCaseStrings;

    // This adds all the valid aliases to aliasesToSwitchStrings.
    void InitializeSlhaOneOrTwoAliases();

    // This adds newParameter to activeDerivedParameters and updates
    // numberOfDistinctActiveParameters and activeParametersToIndices (all
    // aliases in aliasesToCaseStrings which map to caseString are added
    // mapping to the index of this parameter), and returns true paired with
    // the index of newParameter.
    std::pair< bool, size_t >
    AddNewDerivedParameter( std::string const& caseString,
                            SlhaSourcedParameterFunctionoid* newParameter );

    // This checks parameterName against all the special cases. If the
    // parameter name is not recognized as a valid special case, then the
    // LesHouchesAccordBlockEntryManager result is returned, otherwise the
    // special case is registered and returned.
    virtual std::pair< bool, size_t >
    RegisterUnregisteredParameter( std::string const& parameterName );

    // This adds the parameter based on the alias given by switchString for the
    // parameter.
    virtual std::pair< bool, size_t >
    RegisterUnregisteredSpecialCase( std::string const& caseString );

    // This registers the required TX[0,0], AX[0,0], and YX[0,0] parameters
    // where X is replaced by sfermionType and 0 by generationIndex, and then
    // returns the result of AddNewDerivedParameter using a new
    // SlhaTrilinearDiagonalFunctionoid constructed with those parameter
    // functionoids.
    std::pair< bool, size_t >
    RegisterSlhaOneOrTwoCompatibleTrilinear( std::string const& caseString,
                                             char const sfermionType,
                                             char const generationIndex );

    // This registers the required M2X[0,0] and MSOFT[msoftIndex] parameters
    // where X is replaced by sfermionType and 0 by generationIndex, and then
    // returns the result of AddNewDerivedParameter using a new
    // SlhaMassSquaredDiagonalFunctionoid constructed with those parameter
    // functionoids.
    std::pair< bool, size_t >
    RegisterSlhaOneOrTwoCompatibleMassSquared( std::string const& caseString,
                                               char const sfermionType,
                                               char const generationIndex,
                                               std::string const& msoftIndex );

    // This maps both caseString and the format-corrected version of
    // blockString to caseString in aliasesToSwitchStrings.
    void MapCaseStringAndSlhaBlockToCaseString( std::string const& caseString,
                                               std::string const& blockString )
    { aliasesToCaseStrings[ caseString ]
      = aliasesToCaseStrings[ FormatVariable( blockString ) ] = caseString; }

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
  // the index of newParameter.
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::AddNewDerivedParameter(
                                                 std::string const& caseString,
                                SlhaSourcedParameterFunctionoid* newParameter )
  {
    activeDerivedParameters.push_back( newParameter );
    ++numberOfDistinctActiveParameters;
    for( std::map< std::string, std::string >::const_iterator
         aliasToSwitchString( aliasesToCaseStrings.begin() );
         aliasToSwitchString != aliasesToCaseStrings.end();
         ++aliasToSwitchString )
    {
      if( aliasToSwitchString->second == caseString )
      {
        activeParametersToIndices[ aliasToSwitchString->first ]
        = newParameter->IndexInValuesVector();
      }
    }
    return std::pair< bool, size_t >( true,
                                      newParameter->IndexInValuesVector() );
  }

  // This checks parameterName against all the special cases. If the
  // parameter name is not recognized as a valid special case, then the
  // LesHouchesAccordBlockEntryManager result is returned, otherwise the
  // special case is registered and returned.
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::RegisterUnregisteredParameter(
                                             std::string const& parameterName )
  {
    std::map< std::string, std::string >::const_iterator
    aliasToSwitchString( aliasesToCaseStrings.find( parameterName ) );
    if( aliasToSwitchString != aliasesToCaseStrings.end() )
    {
      return RegisterUnregisteredSpecialCase( aliasToSwitchString->second );
    }
    else
    {
      return LesHouchesAccordBlockEntryManager::RegisterUnregisteredParameter(
                                                               parameterName );
    }
  }


  // This registers the required TX[0,0], AX[0,0], and YX[0,0] parameters
  // where X is replaced by sfermionType and 0 by generationIndex, and then
  // returns the result of AddNewDerivedParameter using a new
  // SlhaTrilinearDiagonalFunctionoid constructed with those parameter
  // functionoids.
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::RegisterSlhaOneOrTwoCompatibleTrilinear(
                                                 std::string const& caseString,
                                                       char const sfermionType,
                                                   char const generationIndex )
  {
    std::string blockNameWithIndices( "T?[?,?]" );
    blockNameWithIndices[ 1 ] = sfermionType;
    blockNameWithIndices[ 3 ] = blockNameWithIndices[ 5 ] = generationIndex;
    SlhaSourcedParameterFunctionoid const& directTrilinear(
                RegisterBlockEntry( FormatVariable( blockNameWithIndices ) ) );
    blockNameWithIndices[ 0 ] = 'A';
    SlhaSourcedParameterFunctionoid const& trilinearOverYukawa(
                RegisterBlockEntry( FormatVariable( blockNameWithIndices ) ) );
    blockNameWithIndices[ 0 ] = 'Y';
    SlhaSourcedParameterFunctionoid const& appropriateYukawa(
                RegisterBlockEntry( FormatVariable( blockNameWithIndices ) ) );
    return AddNewDerivedParameter( caseString,
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
                                                 std::string const& caseString,
                                                       char const sfermionType,
                                                    char const generationIndex,
                                                std::string const& msoftIndex )
  {
    std::string blockNameWithIndices( "MS?2[?,?]" );
    blockNameWithIndices[ 2 ] = sfermionType;
    blockNameWithIndices[ 5 ] = blockNameWithIndices[ 7 ] = generationIndex;
    SlhaSourcedParameterFunctionoid const& squareMass(
                RegisterBlockEntry( FormatVariable( blockNameWithIndices ) ) );
    blockNameWithIndices.assign( "MSOFT[" );
    blockNameWithIndices.append( msoftIndex );
    blockNameWithIndices.append( "]" );
    SlhaSourcedParameterFunctionoid const&
    linearMass( RegisterBlockEntry( FormatVariable( blockNameWithIndices ) ) );
    return AddNewDerivedParameter( caseString,
                                   new SlhaMassSquaredDiagonalFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                                    squareMass,
                                                                linearMass ) );
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
