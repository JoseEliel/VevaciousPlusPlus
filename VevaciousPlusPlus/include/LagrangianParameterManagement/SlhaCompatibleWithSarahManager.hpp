/*
 * SlhaCompatibleWithSarahManager.hpp
 *
 *  Created on: Nov 3, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHACOMPATIBLEWITHSARAHMANAGER_HPP_
#define SLHACOMPATIBLEWITHSARAHMANAGER_HPP_

#include "CommonIncludes.hpp"
#include "SlhaBlocksWithSpecialCasesManager.hpp"
#include "SlhaSourcedParameterFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaTwoSourceFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaDsbHiggsVevFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaHiggsMixingBilinearFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaCompatibleWithSarahManager :
                                       public SlhaBlocksWithSpecialCasesManager
  {
  public:
    SlhaCompatibleWithSarahManager( std::string const& validBlocksString,
                                    std::string const& minimumScaleType,
                                    std::string const& minimumScaleArgument,
                                    std::string const& fixedScaleType,
                                    std::string const& fixedScaleArgument );
    SlhaCompatibleWithSarahManager(
                                 std::set< std::string > const& validBlocksSet,
                                    std::string const& minimumScaleType,
                                    std::string const& minimumScaleArgument,
                                    std::string const& fixedScaleType,
                                    std::string const& fixedScaleArgument );
    virtual ~SlhaCompatibleWithSarahManager();


  protected:
    // This adds all the valid aliases to aliasesToSwitchStrings.
    void InitializeSarahAliases()
    {
      // Putting the simple block names as special cases for SARAH means that
      // a SARAH-generated model file which assumes the extra SARAH blocks will
      // still work with non-SARAH SLHA files, to the extent that the DRbar
      // values will be used as both the SARAH tree and loop values.
      CoverAllCasesForSlhaBlock( "DsbVd",
                                 "HMIX[102]" );
      CoverAllCasesForSlhaBlock( "DsbVu",
                                 "HMIX[103]" );
      // The constructor for the base SlhaBlocksWithSpecialCasesManager covers
      // adding the basic Bmu to the alias mapping, even though this derived
      // class over-writes what ends up as its functionoid.
      CoverAllCasesForSlhaBlock( "Bmu",
                                 "HMIX[101]" );
      CoverAllCasesForSlhaBlock( "muTree",
                                 "TREEHMIX[1]" );
      CoverAllCasesForSlhaBlock( "muLoop",
                                 "LOOPHMIX[1]" );
      CoverAllCasesForSlhaBlock( "BmuTree",
                                 "TREEHMIX[101]" );
      CoverAllCasesForSlhaBlock( "BmuLoop",
                                 "LOOPHMIX[101]" );
      CoverAllCasesForSlhaBlock( "mHdSqTree",
                                 "TREEMSOFT[ 21 ]" );
      CoverAllCasesForSlhaBlock( "mHdSqLoop",
                                 "LOOPMSOFT[ 21 ]" );
      CoverAllCasesForSlhaBlock( "mHuSqTree",
                                 "TREEMSOFT[ 22 ]" );
      CoverAllCasesForSlhaBlock( "mHuSqLoop",
                                 "LOOPMSOFT[ 22 ]" );
    }

    // This checks parameterName against all the special cases. If the
    // parameter name is not recognized as a valid special case, then the
    // returned pair is false paired with (size_t)(-1).
    virtual std::pair< bool, size_t >
    RegisterUnregisteredParameter( std::string const& parameterName );

    // This adds the parameter based on the alias given by switchString for the
    // parameter, returning either a found SARAH special case, or otherwise the
    // result of SlhaBlocksWithSpecialCasesManager::SlhaOneOrTwoSpecialCase.
    std::pair< bool, size_t > SarahSpecialCaseDefaultingToSlhaOneOrTwo(
                                               std::string const& caseString );

    // This adds a special case for switchString to be mapped to a
    // SlhaTwoSourceFunctionoid looking preferentially at an interpolated block
    // called (sarahPrefix + slhaName) then defaulting to the block called just
    // slhaName.
    std::pair< bool, size_t >
    RegisterSarahPrefixedSlhaBlock( std::string const& caseString,
                                    std::string const& sarahPrefix,
                                    std::string const& slhaName );
  };





  // This checks parameterName against all the special cases. If the
  // parameter name is not recognized as a valid special case, then the
  // returned pair is false paired with (size_t)(-1).
  inline std::pair< bool, size_t >
  SlhaCompatibleWithSarahManager::RegisterUnregisteredParameter(
                                             std::string const& parameterName )
  {
    std::map< std::string, std::string >::const_iterator
    aliasToSwitchString( aliasesToCaseStrings.find( parameterName ) );
    if( aliasToSwitchString == aliasesToCaseStrings.end() )
    {
      return std::pair< bool, size_t >( false,
                                        -1 );
    }
    return
      SarahSpecialCaseDefaultingToSlhaOneOrTwo( aliasToSwitchString->second );
  }

  // This adds a special case for switchString to be mapped to a
  // SlhaTwoSourceFunctionoid looking preferentially at an interpolated block
  // called (sarahPrefix + slhaName) then defaulting to the block called just
  // slhaName.
  inline std::pair< bool, size_t >
  SlhaCompatibleWithSarahManager::RegisterSarahPrefixedSlhaBlock(
                                                 std::string const& caseString,
                                                std::string const& sarahPrefix,
                                                  std::string const& slhaName )
  {
    SlhaSourcedParameterFunctionoid const&
    sarahFunctionoid( RegisterBlockEntry( sarahPrefix + slhaName ) );
    SlhaSourcedParameterFunctionoid const&
    shouldBeDrbarFunctionoid( RegisterBlockEntry( slhaName ) );
    return AddNewDerivedParameter( caseString,
                new SlhaTwoSourceFunctionoid( numberOfDistinctActiveParameters,
                                              sarahFunctionoid,
                                              shouldBeDrbarFunctionoid ) );
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHACOMPATIBLEWITHSARAHMANAGER_HPP_ */
