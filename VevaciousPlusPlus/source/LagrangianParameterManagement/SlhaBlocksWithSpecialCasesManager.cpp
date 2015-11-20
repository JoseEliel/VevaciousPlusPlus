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
                                         std::string const& fixedScaleArgument,
                                           std::string const& maximumScaleType,
                                    std::string const& maximumScaleArgument ) :
    LesHouchesAccordBlockEntryManager( validBlocksString,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument,
                                       maximumScaleType,
                                       maximumScaleArgument ),
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
                                         std::string const& fixedScaleArgument,
                                           std::string const& maximumScaleType,
                                    std::string const& maximumScaleArgument ) :
    LesHouchesAccordBlockEntryManager( validBlocks,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument,
                                       maximumScaleType,
                                       maximumScaleArgument ),
    activeDerivedParameters(),
    aliasesToCaseStrings()
  {
    InitializeSlhaOneOrTwoAliases();
  }

  SlhaBlocksWithSpecialCasesManager::SlhaBlocksWithSpecialCasesManager(
                                             std::string const& xmlFileName ) :
    LesHouchesAccordBlockEntryManager( xmlFileName ),
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


  // This duplicates a lot of code from RegisterUnregisteredSpecialCase, but
  // there doesn't seem to be an elegant way of using the common code as
  // there is too much entanglement with registering new parameters or not.
  double SlhaBlocksWithSpecialCasesManager::OnceOffSpecialCase(
                                                 std::string const& caseString,
                                          double const logarithmOfScale ) const
  {
    if( ( caseString == "DsbVd" ) || ( caseString == "DsbVu" ) )
    {
      SlhaDsbHiggsVevFunctionoid temporaryParameter( 0,
                                                     0,
                                                     0,
                                                   ( caseString == "DsbVu" ) );
      return temporaryParameter( OnceOffParameter( "HMIX[ 3 ]",
                                                   logarithmOfScale ),
                                 OnceOffParameter( "HMIX[ 2 ]",
                                                   logarithmOfScale ) );
    }
    else if( caseString == "Bmu" )
    {
      SlhaHiggsMixingBilinearFunctionoid temporaryParameter( 0,
                                                             0,
                                                             0 );
      return temporaryParameter( OnceOffParameter( "HMIX[ 4 ]",
                                                   logarithmOfScale ),
                                 OnceOffParameter( "HMIX[ 2 ]",
                                                   logarithmOfScale ) );
    }
    else if( ( caseString == "Te11" )
             ||
             ( caseString == "Te22" )
             ||
             ( caseString == "Te33" )
             ||
             ( caseString == "Td11" )
             ||
             ( caseString == "Td22" )
             ||
             ( caseString == "Td33" )
             ||
             ( caseString == "Tu11" )
             ||
             ( caseString == "Tu22" )
             ||
             ( caseString == "Tu33" )
             ||
             ( caseString == "Te33" )
             ||
             ( caseString == "Te33" )
             ||
             ( caseString == "Te33" )
             ||
             ( caseString == "Te33" )
             ||
             ( caseString == "Te33" ) )
    {
      return OnceOffSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                    toupper( caseString[ 1 ] ),
                                                     caseString[ 2 ],
                                                     logarithmOfScale );
    }
    else if( ( caseString == "Msl211" )
             ||
             ( caseString == "Msl222" )
             ||
             ( caseString == "Msl233" )
             ||
             ( caseString == "Mse211" )
             ||
             ( caseString == "Mse222" )
             ||
             ( caseString == "Mse233" )
             ||
             ( caseString == "Msq211" )
             ||
             ( caseString == "Msq222" )
             ||
             ( caseString == "Msq233" )
             ||
             ( caseString == "Msu211" )
             ||
             ( caseString == "Msu222" )
             ||
             ( caseString == "Msu233" )
             ||
             ( caseString == "Msd211" )
             ||
             ( caseString == "Msd222" )
             ||
             ( caseString == "Msd233" ) )
    {
      // First determine the MSOFT index of the 1st generation, assuming that
      // it is the L doublet and then checking for other sfermion types.
      unsigned int msoftIndex( 31 );
      if( caseString[ 2 ] == 'e' )
      {
        msoftIndex = 34;
      }
      else if( caseString[ 2 ] == 'q' )
      {
        msoftIndex = 41;
      }
      else if( caseString[ 2 ] == 'u' )
      {
        msoftIndex = 44;
      }
      else if( caseString[ 2 ] == 'd' )
      {
        msoftIndex = 47;
      }
      // Next add 1 or 2 if the generation index is '2' or '3' respectively.
      if( caseString[ 4 ] == '2' )
      {
        msoftIndex += 1;
      }
      else if( caseString[ 4 ] == '3' )
      {
        msoftIndex += 2;
      }
      return OnceOffSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                    toupper( caseString[ 2 ] ),
                                                       caseString[ 4 ],
                                                       msoftIndex,
                                                       logarithmOfScale );
    }
    else
    {
      std::stringstream errorBuilder;
      errorBuilder
      << "SlhaBlocksWithSpecialCasesManager::OnceOffSpecialCase could not find"
      << " a case corresponding to \"" << caseString << "\".";
      throw std::runtime_error( errorBuilder.str() );
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
      << "[ \"" << aliasToCase->first << "\", \""
      << aliasToCase->second << "\" ]";
    }
    stringBuilder << " }";

    std::list< SlhaSourcedParameterFunctionoid const* > parameterSortingList;
    typedef LesHouchesAccordBlockEntryManager::LhaBlockEntryInterpolator*
            InterpolatorPointer;
    for( std::vector< InterpolatorPointer >::const_iterator
         baseInterpolator( referenceSafeActiveParameters.begin() );
         baseInterpolator < referenceSafeActiveParameters.end();
         ++baseInterpolator )
    {
      parameterSortingList.push_back( *baseInterpolator );
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
    outputLogScale( log( AppropriateSingleFixedScale() ) );
    std::vector< double > const
    parameterValues( ParameterValues( outputLogScale ) );
    stringBuilder
    << std::endl << "Sorted parameters at scale "
    << AppropriateSingleFixedScale() << " (log = " << outputLogScale
    << ") = {" << std::endl;
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
      << ": \"" << sortedNameAndParameter->first << "\" = "
      << (*sortedNameAndParameter->second)( outputLogScale,
                                            parameterValues );
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
      size_t const vevIndex( RegisterParameter( "HMIX[ 3 ]" ).second );
      size_t const tanBetaIndex( RegisterParameter( "HMIX[ 2 ]" ).second );
      return AddNewDerivedParameter( caseString,
                                     new SlhaDsbHiggsVevFunctionoid(
                                             numberOfDistinctActiveParameters,
                                                                     vevIndex,
                                                                  tanBetaIndex,
                                                 ( caseString == "DsbVu" ) ) );
    }
    else if( caseString == "Bmu" )
    {
      size_t const treePseudoscalarMassSquaredIndex(
                                     RegisterParameter( "HMIX[ 4 ]" ).second );
      size_t const tanBetaIndex( RegisterParameter( "HMIX[ 2 ]" ).second );
      return AddNewDerivedParameter( caseString,
                                     new SlhaHiggsMixingBilinearFunctionoid(
                                              numberOfDistinctActiveParameters,
                                              treePseudoscalarMassSquaredIndex,
                                                              tanBetaIndex ) );
    }
    else if( ( caseString == "Te11" )
             ||
             ( caseString == "Te22" )
             ||
             ( caseString == "Te33" )
             ||
             ( caseString == "Td11" )
             ||
             ( caseString == "Td22" )
             ||
             ( caseString == "Td33" )
             ||
             ( caseString == "Tu11" )
             ||
             ( caseString == "Tu22" )
             ||
             ( caseString == "Tu33" )
             ||
             ( caseString == "Te33" )
             ||
             ( caseString == "Te33" )
             ||
             ( caseString == "Te33" )
             ||
             ( caseString == "Te33" )
             ||
             ( caseString == "Te33" ) )
    {
      return RegisterSlhaOneOrTwoCompatibleTrilinear( caseString,
                                                    toupper( caseString[ 1 ] ),
                                                      caseString[ 2 ] );
    }
    else if( ( caseString == "Msl211" )
             ||
             ( caseString == "Msl222" )
             ||
             ( caseString == "Msl233" )
             ||
             ( caseString == "Mse211" )
             ||
             ( caseString == "Mse222" )
             ||
             ( caseString == "Mse233" )
             ||
             ( caseString == "Msq211" )
             ||
             ( caseString == "Msq222" )
             ||
             ( caseString == "Msq233" )
             ||
             ( caseString == "Msu211" )
             ||
             ( caseString == "Msu222" )
             ||
             ( caseString == "Msu233" )
             ||
             ( caseString == "Msd211" )
             ||
             ( caseString == "Msd222" )
             ||
             ( caseString == "Msd233" ) )
    {
      // First determine the MSOFT index of the 1st generation, assuming that
      // it is the L doublet and then checking for other sfermion types.
      unsigned int msoftIndex( 31 );
      if( caseString[ 2 ] == 'e' )
      {
        msoftIndex = 34;
      }
      else if( caseString[ 2 ] == 'q' )
      {
        msoftIndex = 41;
      }
      else if( caseString[ 2 ] == 'u' )
      {
        msoftIndex = 44;
      }
      else if( caseString[ 2 ] == 'd' )
      {
        msoftIndex = 47;
      }
      // Next add 1 or 2 if the generation index is '2' or '3' respectively.
      if( caseString[ 4 ] == '2' )
      {
        msoftIndex += 1;
      }
      else if( caseString[ 4 ] == '3' )
      {
        msoftIndex += 2;
      }
      return RegisterSlhaOneOrTwoCompatibleMassSquared( caseString,
                                                    toupper( caseString[ 2 ] ),
                                                        caseString[ 4 ],
                                                        msoftIndex );
    }
    else
    {
      return std::pair< bool, size_t >( false,
                                        -1 );
    }
  }

} /* namespace VevaciousPlusPlus */
