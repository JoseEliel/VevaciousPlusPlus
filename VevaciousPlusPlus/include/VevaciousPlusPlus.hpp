/*
 * VevaciousPlusPlus.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef VEVACIOUSPLUSPLUS_HPP_
#define VEVACIOUSPLUSPLUS_HPP_

#include "StandardIncludes.hpp"
#include "PotentialEvaluation.hpp"
#include "PotentialMinimization.hpp"
#include "TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{
  class VevaciousPlusPlus
  {
  public:
    VevaciousPlusPlus( BOL::ArgumentParser& argumentParser,
                       SlhaManager& slhaManager,
                       PotentialMinimizer& potentialMinimizer,
                       TunnelingCalculator& tunnelingCalculator );
    virtual
    ~VevaciousPlusPlus();

    void RunPoint( std::string const& parameterFilename );

    void WriteXmlResults( std::string const& xmlFilename );

    void WriteSlhaResults( std::string const& slhaFilename,
                           bool const writeWarnings = true );


  protected:
    static std::string const versionString;
    static std::string const citationString;

    // BOL::BasicTimer runTimer;
    // It's too much effort to put in a "try to quit as soon as possible after
    // this many seconds per point" functionality - too much chance of a
    // prematurely-ended calculation being misinterpreted by the user as a
    // final result.
    SlhaManager& slhaManager;
    PotentialMinimizer& potentialMinimizer;
    TunnelingCalculator& tunnelingCalculator;
  };

} /* namespace VevaciousPlusPlus */
#endif /* VEVACIOUSPLUSPLUS_HPP_ */
