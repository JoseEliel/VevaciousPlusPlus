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
                                       std::string const& validBlocksString ) :
    SlhaBlocksWithSpecialCasesManager( validBlocksString )
  {
    InitializeSarahAliases();
  }

  SlhaCompatibleWithSarahManager::SlhaCompatibleWithSarahManager(
                                 std::set< std::string > const& validBlocks ) :
    SlhaBlocksWithSpecialCasesManager( validBlocks )
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
                                              std::string const& switchString )
  {
    switch( switchString )
    {
      case "DsbVd":
      case "DsbVu":
        bool upNotDown( switchString == "DsbVu" );
        SlhaSourcedParameterFunctionoid const& sarahFunctionoid(
             RegisterBlockEntry( upNotDown ? "HMIX[ 103 ]" : "HMIX[ 102 ]" ) );
        std::pair< bool, size_t >
        slhaResult( SlhaOneOrTwoSpecialCase( switchString ) );
        SlhaSourcedParameterFunctionoid const&
        slhaFunctionoid( *(FindActiveDerivedParameter( slhaResult.second )) );
        // This will overwrite activeParametersToIndices[ switchString ] to
        // map switchString to numberOfDistinctActiveParameters as its index in
        // the values vector (corresponding now to the new
        // SlhaTwoSourceFunctionoid), rather than to slhaFunctionoid as was set
        // by SlhaOneOrTwoSpecialCase( switchString ).
        return AddNewDerivedParameter( switchString,
                new SlhaTwoSourceFunctionoid( numberOfDistinctActiveParameters,
                                              sarahFunctionoid,
                                              slhaFunctionoid ) );
        break;
      case "Bmu":
        SlhaSourcedParameterFunctionoid const&
        sarahFunctionoid( RegisterBlockEntry( "HMIX[ 101 ]" ) );
        std::pair< bool, size_t >
        slhaResult( SlhaOneOrTwoSpecialCase( switchString ) );
        SlhaSourcedParameterFunctionoid const&
        slhaFunctionoid( *(FindActiveDerivedParameter( slhaResult.second )) );
        // This will overwrite activeParametersToIndices[ switchString ] to
        // map switchString to numberOfDistinctActiveParameters as its index in
        // the values vector (corresponding now to the new
        // SlhaTwoSourceFunctionoid), rather than to slhaFunctionoid as was set
        // by SlhaOneOrTwoSpecialCase( switchString ).
        return AddNewDerivedParameter( switchString,
                new SlhaTwoSourceFunctionoid( numberOfDistinctActiveParameters,
                                              sarahFunctionoid,
                                              slhaFunctionoid ) );
        break;
      case "muTree":
        return RegisterSarahPrefixedSlhaBlock( switchString,
                                               "TREE",
                                               "HMIX[ 1 ]" );
        break;
      case "muLoop":
        return RegisterSarahPrefixedSlhaBlock( switchString,
                                               "LOOP",
                                               "HMIX[ 1 ]" );
        break;
      case "BmuTree":
        SlhaSourcedParameterFunctionoid const&
        sarahFunctionoid( RegisterBlockEntry( "TREEHMIX[ 101 ]" ) );
        std::pair< bool, size_t > shouldBeDrbarResult(
                           SarahSpecialCaseDefaultingToSlhaOneOrTwo( "Bmu" ) );
        SlhaSourcedParameterFunctionoid const& shouldBeDrbarFunctionoid(
                 *(FindActiveDerivedParameter( shouldBeDrbarResult.second )) );
        // This will overwrite activeParametersToIndices[ switchString ] to
        // map switchString to numberOfDistinctActiveParameters as its index in
        // the values vector (corresponding now to the new
        // SlhaTwoSourceFunctionoid), rather than to shouldBeDrbarFunctionoid
        // as was set by SarahSpecialCaseDefaultingToSlhaOneOrTwo( "Bmu" ).
        return AddNewDerivedParameter( switchString,
                new SlhaTwoSourceFunctionoid( numberOfDistinctActiveParameters,
                                              sarahFunctionoid,
                                              shouldBeDrbarFunctionoid ) );
        break;
      case "BmuLoop":
        SlhaSourcedParameterFunctionoid const&
        sarahFunctionoid( RegisterBlockEntry( "LOOPHMIX[ 101 ]" ) );
        std::pair< bool, size_t > shouldBeDrbarResult(
                           SarahSpecialCaseDefaultingToSlhaOneOrTwo( "Bmu" ) );
        SlhaSourcedParameterFunctionoid const& shouldBeDrbarFunctionoid(
                 *(FindActiveDerivedParameter( shouldBeDrbarResult.second )) );
        // This will overwrite activeParametersToIndices[ switchString ] to
        // map switchString to numberOfDistinctActiveParameters as its index in
        // the values vector (corresponding now to the new
        // SlhaTwoSourceFunctionoid), rather than to shouldBeDrbarFunctionoid
        // as was set by SarahSpecialCaseDefaultingToSlhaOneOrTwo( "Bmu" ).
        return AddNewDerivedParameter( switchString,
                new SlhaTwoSourceFunctionoid( numberOfDistinctActiveParameters,
                                              sarahFunctionoid,
                                              shouldBeDrbarFunctionoid ) );
        break;
      case "mHdSqTree":
        return RegisterSarahPrefixedSlhaBlock( switchString,
                                               "TREE",
                                               "MSOFT[ 21 ]" );
        break;
      case "mHdSqLoop":
        return RegisterSarahPrefixedSlhaBlock( switchString,
                                               "LOOP",
                                               "MSOFT[ 21 ]" );
        break;
      case "mHuSqTree":
        return RegisterSarahPrefixedSlhaBlock( switchString,
                                               "TREE",
                                               "MSOFT[ 22 ]" );
        break;
      case "mHuSqLoop":
        return RegisterSarahPrefixedSlhaBlock( switchString,
                                               "LOOP",
                                               "MSOFT[ 22 ]" );
        break;
      default:
        return SlhaOneOrTwoSpecialCase( switchString );
    }
  }

} /* namespace VevaciousPlusPlus */
