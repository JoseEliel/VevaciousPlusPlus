/*
 * SlhaBlocksWithSpecialCasesManager.cpp
 *
 *  Created on: Nov 2, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaBlocksWithSpecialCasesManager.hpp"

namespace VevaciousPlusPlus
{

  SlhaBlocksWithSpecialCasesManager::SlhaBlocksWithSpecialCasesManager(
                                          std::string const& validBlocksString,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                      std::string const& fixedScaleArgument ) :
    LesHouchesAccordBlockEntryManager( validBlocksString,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument ),
    activeDerivedParameters(),
    aliasesToCaseStrings()
  {
    InitializeSlhaOneOrTwoAliases();
  }

  SlhaBlocksWithSpecialCasesManager::SlhaBlocksWithSpecialCasesManager(
                                 std::set< std::string > const& validBlocksSet,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                      std::string const& fixedScaleArgument ) :
    LesHouchesAccordBlockEntryManager( validBlocks,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument ),
    activeDerivedParameters(),
    aliasesToCaseStrings()
  {
    InitializeSlhaOneOrTwoAliases();
  }

  SlhaBlocksWithSpecialCasesManager::~SlhaBlocksWithSpecialCasesManager()
  {
    for( size_t deletionIndex( 0 );
         deletionIndex < activeDerivedParameters.size();
         ++deletionIndex )
    {
      delete activeDerivedParameters[ deletionIndex ];
    }
  }


  // This is mainly for debugging.
  std::string SlhaBlocksWithSpecialCasesManager::AsDebuggingString() const
  {
    std::stringstream stringBuilder;
    stringBuilder << LesHouchesAccordBlockEntryManager::AsDebuggingString()
    << std::endl;
    stringBuilder << "activeDerivedParameters = { " << std::endl;
    for( std::vector< SlhaSourcedParameterFunctionoid* >::const_iterator
         derivedParameter( activeDerivedParameters.begin() );
         derivedParameter < activeDerivedParameters.end();
         ++derivedParameter )
    {
      stringBuilder << (*derivedParameter)->AsDebuggingString() << std::endl;
    }
    stringBuilder << "}" << std::endl << "aliasesToCaseStrings = { ";
    for( std::map< std::string, std::string >::const_iterator
         aliasToCase( aliasesToCaseStrings.begin() );
         aliasToCase != aliasesToCaseStrings.end();
         ++aliasToCase )
    {
      if( aliasToCase != aliasesToCaseStrings.begin() )
      {
        stringBuilder << "," << std::endl;
      }
      stringBuilder
      << "[ \"" << aliasToCase->first << "\", "
      << aliasToCase->second << " ]";
    }
    stringBuilder << " }";

    std::list< SlhaSourcedParameterFunctionoid const* > parameterSortingList;
    typedef LesHouchesAccordBlockEntryManager::LhaBlockEntryInterpolator
            BaseInterpolator;
    for( std::vector< BaseInterpolator >::const_iterator
         baseInterpolator( referenceSafeActiveParameters.begin() );
         baseInterpolator < referenceSafeActiveParameters.end();
         ++baseInterpolator )
    {
      parameterSortingList.push_back( &(*baseInterpolator) );
    }
    for( std::vector< SlhaSourcedParameterFunctionoid* >::const_iterator
         derivedParameter( activeDerivedParameters.begin() );
         derivedParameter < activeDerivedParameters.end();
         ++derivedParameter )
    {
      parameterSortingList.push_back( *derivedParameter );
    }
    parameterSortingList.sort( &SortParameterByIndex );
    typedef std::pair< std::string, SlhaSourcedParameterFunctionoid const* >
    nameAndParameter;
    std::vector< nameAndParameter > sortedNamesAndParameters;
    for( std::list< SlhaSourcedParameterFunctionoid const* >::const_iterator
         sortedParameter( parameterSortingList.begin() );
         sortedParameter != parameterSortingList.end();
         ++sortedParameter )
    {
      sortedNamesAndParameters.push_back( std::make_pair( "unknown",
                                                          *sortedParameter ) );
    }
    for( std::map< std::string, size_t >::const_iterator
         parameterToIndex( activeParametersToIndices.begin() );
         parameterToIndex != activeParametersToIndices.end();
         ++parameterToIndex )
    {
      std::string& parameterName(
                  sortedNamesAndParameters[ parameterToIndex->second ].first );
      if( parameterName == "unknown" )
      {
        parameterName = parameterToIndex->first;
      }
      else
      {
        parameterName.append( ", " ).append( parameterToIndex->first );
      }
    }
    double const
    outputLogScale( log( AppropriateFixedScaleForParameterPoint() ) );
    stringBuilder
    << std::endl << "Sorted parameters at scale "
    << AppropriateFixedScaleForParameterPoint() << " (log = " << outputLogScale
    << ")= {" << std::endl;
    for( std::vector< nameAndParameter >::const_iterator
         sortedNameAndParameter( sortedNamesAndParameters.begin() );
         sortedNameAndParameter < sortedNamesAndParameters.end();
         ++sortedNameAndParameter )
    {
      if( sortedNameAndParameter != sortedNamesAndParameters.begin() )
      {
        stringBuilder << "," << std::endl;
      }

      stringBuilder  << sortedNameAndParameter->second->IndexInValuesVector()
      << ": \"" << sortedNameAndParameter->first << "\" = ";
      // debugging:
      /**/
      if( sortedNameAndParameter->second->IndexInValuesVector() < 15 )
      {
        stringBuilder << (*sortedNameAndParameter->second)( outputLogScale );
      }
      else
      {
        stringBuilder << "not available until I sort out these bugs";
      }
      /**/
    }
    stringBuilder << " }" << std::endl;
    return stringBuilder.str();
  }

  // This adds all the valid aliases to aliasesToSwitchStrings.
  void SlhaBlocksWithSpecialCasesManager::InitializeSlhaOneOrTwoAliases()
  {
    aliasesToCaseStrings[ "DsbVd" ] = "DsbVd";
    aliasesToCaseStrings[ "DsbVu" ] = "DsbVu";
    aliasesToCaseStrings[ "Bmu" ]
    = aliasesToCaseStrings[ "m3Sq" ] = "Bmu";
    MapCaseStringAndSlhaBlockToCaseString( "Te11",
                                           "TE[1,1]" );
    MapCaseStringAndSlhaBlockToCaseString( "Te22",
                                           "TE[2,2]" );
    MapCaseStringAndSlhaBlockToCaseString( "Te33",
                                           "TE[3,3]" );
    MapCaseStringAndSlhaBlockToCaseString( "Td11",
                                           "TD[1,1]" );
    MapCaseStringAndSlhaBlockToCaseString( "Td22",
                                           "TD[2,2]" );
    MapCaseStringAndSlhaBlockToCaseString( "Td33",
                                           "TD[3,3]" );
    MapCaseStringAndSlhaBlockToCaseString( "Tu11",
                                           "TU[1,1]" );
    MapCaseStringAndSlhaBlockToCaseString( "Tu22",
                                           "TU[2,2]" );
    MapCaseStringAndSlhaBlockToCaseString( "Tu33",
                                           "TU[3,3]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msl211",
                                           "MSL2[1,1]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msl222",
                                           "MSL2[2,2]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msl233",
                                           "MSL2[3,3]" );
    MapCaseStringAndSlhaBlockToCaseString( "Mse211",
                                           "MSE2[1,1]" );
    MapCaseStringAndSlhaBlockToCaseString( "Mse222",
                                           "MSE2[2,2]" );
    MapCaseStringAndSlhaBlockToCaseString( "Mse233",
                                           "MSE2[3,3]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msq211",
                                           "MSQ2[1,1]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msq222",
                                           "MSQ2[2,2]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msq233",
                                           "MSQ2[3,3]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msu211",
                                           "MSU2[1,1]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msu222",
                                           "MSU2[2,2]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msu233",
                                           "MSU2[3,3]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msd211",
                                           "MSD2[1,1]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msd222",
                                           "MSD2[2,2]" );
    MapCaseStringAndSlhaBlockToCaseString( "Msd233",
                                           "MSD2[3,3]" );
  }

  // This adds the parameter based on the alias given by switchString for the
  // parameter.
  std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::RegisterUnregisteredSpecialCase(
                                                std::string const& caseString )
  {
    if( ( caseString == "DsbVd" ) || ( caseString == "DsbVu" ) )
    {
      SlhaSourcedParameterFunctionoid const&
      vevLength( RegisterBlockEntry( FormatVariable( "HMIX[ 3 ]" ) ) );
      SlhaSourcedParameterFunctionoid const&
      tanBeta( RegisterBlockEntry( FormatVariable( "HMIX[ 2 ]" ) ) );

      return AddNewDerivedParameter( caseString,
                                     new SlhaDsbHiggsVevFunctionoid(
                                             numberOfDistinctActiveParameters,
                                                                     vevLength,
                                                                     tanBeta,
                                                 ( caseString == "DsbVu" ) ) );
    }
    else if( caseString == "Bmu" )
    {
      SlhaSourcedParameterFunctionoid const& treePseudoscalarMassSquared(
                         RegisterBlockEntry( FormatVariable( "HMIX[ 4 ]" ) ) );
      SlhaSourcedParameterFunctionoid const&
      tanBeta( RegisterBlockEntry( FormatVariable( "HMIX[ 2 ]" ) ) );
      return AddNewDerivedParameter( caseString,
                                     new SlhaHiggsMixingBilinearFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                   treePseudoscalarMassSquared,
                                                                   tanBeta ) );
    }
    else if( caseString == "Te11" )
    {
      return RegisterSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                      'E',
                                                      '1' );
    }
    else if( caseString == "Te22" )
    {
      return RegisterSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                      'E',
                                                      '2' );
    }
    else if( caseString == "Te33" )
    {
      return RegisterSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                      'E',
                                                      '3' );
    }
    else if( caseString == "Td11" )
    {
      return RegisterSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                      'D',
                                                      '1' );
    }
    else if( caseString == "Td22" )
    {
      return RegisterSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                      'D',
                                                      '2' );
    }
    else if( caseString == "Td33" )
    {
      return RegisterSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                      'D',
                                                      '3' );
    }
    else if( caseString == "Tu11" )
    {
      return RegisterSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                      'U',
                                                      '1' );
    }
    else if( caseString == "Tu22" )
    {
      return RegisterSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                      'U',
                                                      '2' );
    }
    else if( caseString == "Tu33" )
    {
      return RegisterSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                      'U',
                                                      '3' );
    }
    else if( caseString == "Msl211" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'L',
                                                        '1',
                                                        "31" );
    }
    else if( caseString == "Msl222" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'L',
                                                        '2',
                                                        "32" );
    }
    else if( caseString == "Msl233" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'L',
                                                        '3',
                                                        "33" );
    }
    else if( caseString == "Mse211" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'E',
                                                        '1',
                                                        "34" );
    }
    else if( caseString == "Mse222" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'E',
                                                        '2',
                                                        "35" );
    }
    else if( caseString == "Mse233" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'E',
                                                        '3',
                                                        "36" );
    }
    else if( caseString == "Msq211" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                          'Q',
                                                          '1',
                                                          "41" );
    }
    else if( caseString == "Msq222" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'Q',
                                                        '2',
                                                        "42" );
    }
    else if( caseString == "Msq233" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'Q',
                                                        '3',
                                                        "43" );
    }
    else if( caseString == "Msu211" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'U',
                                                        '1',
                                                        "44" );
    }
    else if( caseString == "Msu222" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'U',
                                                        '2',
                                                        "45" );
    }
    else if( caseString == "Msu233" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'U',
                                                        '3',
                                                        "46" );
    }
    else if( caseString == "Msd211" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'D',
                                                        '1',
                                                        "47" );
    }
    else if( caseString == "Msd222" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'D',
                                                        '2',
                                                        "48" );
    }
    else if( caseString == "Msd233" )
    {
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                        'D',
                                                        '3',
                                                        "49" );
    }
    else
    {
      return std::pair< bool, size_t >( false,
                                        -1 );
    }
  }

} /* namespace VevaciousPlusPlus */
