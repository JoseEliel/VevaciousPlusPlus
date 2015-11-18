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

    // This duplicates a lot of code from RegisterUnregisteredSpecialCase, but
    // there doesn't seem to be an elegant way of using the common code as
    // there is too much entanglement with registering new parameters or not.
    virtual double OnceOffSpecialCase( std::string const& parameterName,
                                       double const logarithmOfScale ) const;

    // This adds the parameter based on the alias given by switchString for the
    // parameter, returning either a found SARAH special case, or otherwise the
    // result of SlhaBlocksWithSpecialCasesManager::SlhaOneOrTwoSpecialCase.
    virtual std::pair< bool, size_t >
    RegisterUnregisteredSpecialCase(  std::string const& caseString );
  };

} /* namespace VevaciousPlusPlus */

#endif /* SLHACOMPATIBLEWITHSARAHMANAGER_HPP_ */
