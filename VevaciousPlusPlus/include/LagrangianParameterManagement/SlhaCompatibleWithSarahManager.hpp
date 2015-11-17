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
    void InitializeSarahAliases();

    // This adds the parameter based on the alias given by switchString for the
    // parameter, returning either a found SARAH special case, or otherwise the
    // result of SlhaBlocksWithSpecialCasesManager::SlhaOneOrTwoSpecialCase.
    virtual std::pair< bool, size_t >
    RegisterUnregisteredSpecialCase(  std::string const& caseString );

    // This adds a special case for switchString to be mapped to a
    // SlhaTwoSourceFunctionoid looking preferentially at an interpolated block
    // called (sarahPrefix + slhaName) then defaulting to the block called just
    // slhaName.
    std::pair< bool, size_t >
    RegisterSarahPrefixedSlhaBlock( std::string const& caseString,
                                    std::string const& sarahPrefix,
                                    std::string const& slhaName );
  };





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
    SlhaSourcedParameterFunctionoid const& sarahFunctionoid(
              RegisterBlockEntry( FormatVariable( sarahPrefix + slhaName ) ) );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "sarahFunctionoid.AsDebuggingString() =" << std::endl;
    std::cout << sarahFunctionoid.AsDebuggingString() << std::endl;
    std::cout << std::endl;/**/

    SlhaSourcedParameterFunctionoid const& shouldBeDrbarFunctionoid(
                            RegisterBlockEntry( FormatVariable( slhaName ) ) );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "sarahFunctionoid.AsDebuggingString() =" << std::endl;
    std::cout << sarahFunctionoid.AsDebuggingString() << std::endl;
    std::cout << "shouldBeDrbarFunctionoid.AsDebuggingString() =" << std::endl;
    std::cout << shouldBeDrbarFunctionoid.AsDebuggingString() << std::endl;
    std::cout << std::endl;/**/
    return AddNewDerivedParameter( caseString,
                new SlhaTwoSourceFunctionoid( numberOfDistinctActiveParameters,
                                              sarahFunctionoid,
                                              shouldBeDrbarFunctionoid ) );
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHACOMPATIBLEWITHSARAHMANAGER_HPP_ */
