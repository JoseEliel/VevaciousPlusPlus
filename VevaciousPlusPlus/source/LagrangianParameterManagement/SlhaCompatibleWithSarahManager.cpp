/*
 * SlhaCompatibleWithSarahManager.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaCompatibleWithSarahManager.hpp"

namespace VevaciousPlusPlus
{

  SlhaCompatibleWithSarahManager::SlhaCompatibleWithSarahManager(
                                          std::string const& validBlocksString,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                      std::string const& fixedScaleArgument ) :
    SlhaBlocksWithSpecialCasesManager( validBlocksString,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument )
  {
    InitializeSarahAliases();
  }

  SlhaCompatibleWithSarahManager::SlhaCompatibleWithSarahManager(
                                 std::set< std::string > const& validBlocksSet,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                      std::string const& fixedScaleArgument ) :
    SlhaBlocksWithSpecialCasesManager( validBlocks,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument )
  {
    InitializeSarahAliases();
  }

  SlhaCompatibleWithSarahManager::~SlhaCompatibleWithSarahManager()
  {
    // This does nothing beyond what the base SlhaBlocksWithSpecialCasesManager
    // destructor does.
  }


  // This adds the parameter based on the alias given by switchString for the
  // parameter.
  std::pair< bool, size_t >
  SlhaCompatibleWithSarahManager::SarahSpecialCaseDefaultingToSlhaOneOrTwo(
                                                std::string const& caseString )
  {

    if( ( caseString == "DsbVd" ) || ( caseString == "DsbVu" ) )
    {
      bool upNotDown( caseString == "DsbVu" );
      SlhaSourcedParameterFunctionoid const& sarahFunctionoid(
             RegisterBlockEntry( upNotDown ? "HMIX[ 103 ]" : "HMIX[ 102 ]" ) );
      std::pair< bool, size_t >
      slhaResult( SlhaOneOrTwoSpecialCase( caseString ) );
      SlhaSourcedParameterFunctionoid const&
      slhaFunctionoid( *(FindActiveDerivedParameter( slhaResult.second )) );
      // This will overwrite activeParametersToIndices[ switchString ] to
      // map switchString to numberOfDistinctActiveParameters as its index in
      // the values vector (corresponding now to the new
      // SlhaTwoSourceFunctionoid), rather than to slhaFunctionoid as was set
      // by SlhaOneOrTwoSpecialCase( switchString ).
      return AddNewDerivedParameter( caseString,
                new SlhaTwoSourceFunctionoid( numberOfDistinctActiveParameters,
                                              sarahFunctionoid,
                                              slhaFunctionoid ) );
    }
    else if( caseString =="Bmu" )
    {
      SlhaSourcedParameterFunctionoid const&
      sarahFunctionoid( RegisterBlockEntry( "HMIX[ 101 ]" ) );
      std::pair< bool, size_t >
      slhaResult( SlhaOneOrTwoSpecialCase( caseString ) );
      SlhaSourcedParameterFunctionoid const&
      slhaFunctionoid( *(FindActiveDerivedParameter( slhaResult.second )) );
      // This will overwrite activeParametersToIndices[ switchString ] to
      // map switchString to numberOfDistinctActiveParameters as its index in
      // the values vector (corresponding now to the new
      // SlhaTwoSourceFunctionoid), rather than to slhaFunctionoid as was set
      // by SlhaOneOrTwoSpecialCase( switchString ).
      return AddNewDerivedParameter( caseString,
                new SlhaTwoSourceFunctionoid( numberOfDistinctActiveParameters,
                                              sarahFunctionoid,
                                              slhaFunctionoid ) );
    }
    else if( caseString == "muTree" )
    {
      return RegisterSarahPrefixedSlhaBlock( caseString,
                                             "TREE",
                                             "HMIX[ 1 ]" );
    }
    else if( caseString == "muLoop" )
    {
      return RegisterSarahPrefixedSlhaBlock( caseString,
                                             "LOOP",
                                             "HMIX[ 1 ]" );
    }
    else if( caseString == "BmuTree" )
    {
      SlhaSourcedParameterFunctionoid const&
      sarahFunctionoid( RegisterBlockEntry( "TREEHMIX[ 101 ]" ) );
      std::pair< bool, size_t >
      shouldBeDrbarResult( SarahSpecialCaseDefaultingToSlhaOneOrTwo( "Bmu" ) );
      SlhaSourcedParameterFunctionoid const& shouldBeDrbarFunctionoid(
                 *(FindActiveDerivedParameter( shouldBeDrbarResult.second )) );
      // This will overwrite activeParametersToIndices[ switchString ] to
      // map switchString to numberOfDistinctActiveParameters as its index in
      // the values vector (corresponding now to the new
      // SlhaTwoSourceFunctionoid), rather than to shouldBeDrbarFunctionoid
      // as was set by SarahSpecialCaseDefaultingToSlhaOneOrTwo( "Bmu" ).
      return AddNewDerivedParameter( caseString,
                new SlhaTwoSourceFunctionoid( numberOfDistinctActiveParameters,
                                              sarahFunctionoid,
                                              shouldBeDrbarFunctionoid ) );
    }
    else if( caseString == "BmuLoop" )
    {
      SlhaSourcedParameterFunctionoid const&
      sarahFunctionoid( RegisterBlockEntry( "LOOPHMIX[ 101 ]" ) );
      std::pair< bool, size_t >
      shouldBeDrbarResult( SarahSpecialCaseDefaultingToSlhaOneOrTwo( "Bmu" ) );
      SlhaSourcedParameterFunctionoid const& shouldBeDrbarFunctionoid(
                 *(FindActiveDerivedParameter( shouldBeDrbarResult.second )) );
      // This will overwrite activeParametersToIndices[ switchString ] to
      // map switchString to numberOfDistinctActiveParameters as its index in
      // the values vector (corresponding now to the new
      // SlhaTwoSourceFunctionoid), rather than to shouldBeDrbarFunctionoid
      // as was set by SarahSpecialCaseDefaultingToSlhaOneOrTwo( "Bmu" ).
      return AddNewDerivedParameter( caseString,
                new SlhaTwoSourceFunctionoid( numberOfDistinctActiveParameters,
                                              sarahFunctionoid,
                                              shouldBeDrbarFunctionoid ) );
    }
    else if( caseString == "mHdSqTree" )
    {
      return RegisterSarahPrefixedSlhaBlock( caseString,
                                             "TREE",
                                             "MSOFT[ 21 ]" );
    }
    else if( caseString == "mHdSqLoop" )
    {
      return RegisterSarahPrefixedSlhaBlock( caseString,
                                             "LOOP",
                                             "MSOFT[ 21 ]" );
    }
    else if( caseString == "mHuSqTree" )
    {
      return RegisterSarahPrefixedSlhaBlock( caseString,
                                             "TREE",
                                             "MSOFT[ 22 ]" );
    }
    else if( caseString == "mHuSqLoop" )
    {
      return RegisterSarahPrefixedSlhaBlock( caseString,
                                             "LOOP",
                                             "MSOFT[ 22 ]" );
    }
    else
    {
      return SlhaOneOrTwoSpecialCase( caseString );
    }
  }

} /* namespace VevaciousPlusPlus */
