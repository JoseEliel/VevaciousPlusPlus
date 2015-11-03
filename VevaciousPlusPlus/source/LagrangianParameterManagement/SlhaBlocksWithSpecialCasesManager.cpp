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
                                       std::string const& validBlocksString ) :
    LesHouchesAccordBlockEntryManager( validBlocksString ),
    activeDerivedParameters(),
    aliasesToSwitchStrings()
  {
    InitializeSlhaOneOrTwoAliases();
  }

  SlhaBlocksWithSpecialCasesManager::SlhaBlocksWithSpecialCasesManager(
                                 std::set< std::string > const& validBlocks ) :
    LesHouchesAccordBlockEntryManager( validBlocks ),
    activeDerivedParameters(),
    aliasesToSwitchStrings()
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


  // This adds all the valid aliases to aliasesToSwitchStrings.
  void SlhaBlocksWithSpecialCasesManager::InitializeSlhaOneOrTwoAliases()
  {
    aliasesToSwitchStrings[ "DsbVd" ] = "DsbVd";
    aliasesToSwitchStrings[ "DsbVu" ] = "DsbVu";
    aliasesToSwitchStrings[ "Bmu" ]
    = aliasesToSwitchStrings[ "m3Sq" ] = "Bmu";
    CoverAllCasesForSlhaBlock( "Te11",
                               "TE[1,1]" );
    CoverAllCasesForSlhaBlock( "Te22",
                               "TE[2,2]" );
    CoverAllCasesForSlhaBlock( "Te33",
                               "TE[3,3]" );
    CoverAllCasesForSlhaBlock( "Td11",
                               "TD[1,1]" );
    CoverAllCasesForSlhaBlock( "Td22",
                               "TD[2,2]" );
    CoverAllCasesForSlhaBlock( "Td33",
                               "TD[3,3]" );
    CoverAllCasesForSlhaBlock( "Tu11",
                               "TU[1,1]" );
    CoverAllCasesForSlhaBlock( "Tu22",
                               "TU[2,2]" );
    CoverAllCasesForSlhaBlock( "Tu33",
                               "TU[3,3]" );
    CoverAllCasesForSlhaBlock( "Msl211",
                               "MSL2[1,1]" );
    CoverAllCasesForSlhaBlock( "Msl222",
                               "MSL2[2,2]" );
    CoverAllCasesForSlhaBlock( "Msl233",
                               "MSL2[3,3]" );
    CoverAllCasesForSlhaBlock( "Mse211",
                               "MSE2[1,1]" );
    CoverAllCasesForSlhaBlock( "Mse222",
                               "MSE2[2,2]" );
    CoverAllCasesForSlhaBlock( "Mse233",
                               "MSE2[3,3]" );
    CoverAllCasesForSlhaBlock( "Msq211",
                               "MSQ2[1,1]" );
    CoverAllCasesForSlhaBlock( "Msq222",
                               "MSQ2[2,2]" );
    CoverAllCasesForSlhaBlock( "Msq233",
                               "MSQ2[3,3]" );
    CoverAllCasesForSlhaBlock( "Msu211",
                               "MSU2[1,1]" );
    CoverAllCasesForSlhaBlock( "Msu222",
                               "MSU2[2,2]" );
    CoverAllCasesForSlhaBlock( "Msu233",
                               "MSU2[3,3]" );
    CoverAllCasesForSlhaBlock( "Msd211",
                               "MSD2[1,1]" );
    CoverAllCasesForSlhaBlock( "Msd222",
                               "MSD2[2,2]" );
    CoverAllCasesForSlhaBlock( "Msd233",
                               "MSD2[3,3]" );
  }

  // This adds the parameter based on the alias given by switchString for the
  // parameter.
  std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::SlhaOneOrTwoSpecialCase(
                                              std::string const& switchString )
  {
    switch( switchString )
    {
      case "DsbVd":
      case "DsbVu":
        SlhaSourcedParameterFunctionoid const&
        vevLength( RegisterBlockEntry( "HMIX[ 3 ]" ) );
        SlhaSourcedParameterFunctionoid const&
        tanBeta( RegisterBlockEntry( "HMIX[ 2 ]" ) );

        return AddNewDerivedParameter( switchString,
            new SlhaDsbHiggsVevFunctionoid( numberOfDistinctActiveParameters,
                                            vevLength,
                                            tanBeta,
                                            ( switchString == "DsbVu" ) ) );
        break;
      case "Bmu":
        SlhaSourcedParameterFunctionoid const&
        treePseudoscalarMassSquared( RegisterBlockEntry( "HMIX[ 4 ]" ) );
        SlhaSourcedParameterFunctionoid const&
        tanBeta( RegisterBlockEntry( "HMIX[ 2 ]" ) );
        return AddNewDerivedParameter( switchString,
                                      new SlhaHiggsMixingBilinearFunctionoid(
                                            numberOfDistinctActiveParameters,
                                                 treePseudoscalarMassSquared,
                                                                 tanBeta ) );
      case "Te11":
        return RegisterSlhaOneOrTwoCompatibleTrilinear( switchString,
                                                        'E',
                                                        '1' );
        break;
      case "Te22":
        return RegisterSlhaOneOrTwoCompatibleTrilinear( switchString,
                                                        'E',
                                                        '2' );
        break;
      case "Te33":
        return RegisterSlhaOneOrTwoCompatibleTrilinear( switchString,
                                                        'E',
                                                        '3' );
        break;
      case "Td11":
        return RegisterSlhaOneOrTwoCompatibleTrilinear( switchString,
                                                        'D',
                                                        '1' );
        break;
      case "Td22":
        return RegisterSlhaOneOrTwoCompatibleTrilinear( switchString,
                                                        'D',
                                                        '2' );
        break;
      case "Td33":
        return RegisterSlhaOneOrTwoCompatibleTrilinear( switchString,
                                                        'D',
                                                        '3' );
        break;
      case "Tu11":
        return RegisterSlhaOneOrTwoCompatibleTrilinear( switchString,
                                                        'U',
                                                        '1' );
        break;
      case "Tu22":
        return RegisterSlhaOneOrTwoCompatibleTrilinear( switchString,
                                                        'U',
                                                        '2' );
        break;
      case "Tu33":
        return RegisterSlhaOneOrTwoCompatibleTrilinear( switchString,
                                                        'U',
                                                        '3' );
        break;
      case "Msl211":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'L',
                                                         '1',
                                                         "31" );
        break;
      case "Msl222":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'L',
                                                         '1',
                                                         "32" );
        break;
      case "Msl233":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'L',
                                                         '1',
                                                         "33" );
        break;
      case "Mse211":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'E',
                                                         '1',
                                                         "34" );
        break;
      case "Mse222":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'E',
                                                         '2',
                                                         "35" );
        break;
      case "Mse233":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'E',
                                                         '3',
                                                         "36" );
        break;
      case "Msq211":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'Q',
                                                         '1',
                                                         "41" );
        break;
      case "Msq222":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'Q',
                                                         '1',
                                                         "42" );
        break;
      case "Msq233":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'Q',
                                                         '1',
                                                         "43" );
        break;
      case "Msu211":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'U',
                                                         '1',
                                                         "44" );
        break;
      case "Msu222":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'U',
                                                         '2',
                                                         "45" );
        break;
      case "Msu233":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'U',
                                                         '3',
                                                         "46" );
      case "Msd211":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'D',
                                                         '1',
                                                         "47" );
        break;
      case "Msd222":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'D',
                                                         '2',
                                                         "48" );
        break;
      case "Msd233":
        return RegisterSlhaOneOrTwoCompatibleMassSquared( switchString,
                                                         'D',
                                                         '3',
                                                         "49" );
        break;
      default:
        return std::pair< bool, size_t >( false,
                                          -1 );
    }
  }

} /* namespace VevaciousPlusPlus */
