/*
 * VevaciousPlusPlus.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef VEVACIOUSPLUSPLUS_HPP_
#define VEVACIOUSPLUSPLUS_HPP_

#include "StandardIncludes.hpp"
#include "PotentialEvaluation/PotentialEvaluation.hpp"
#include "PotentialMinimization/PotentialMinimization.hpp"
#include "TunnelingCalculation/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{
  class VevaciousPlusPlus
  {
  public:
    VevaciousPlusPlus( BOL::ArgumentParser& argumentParser,
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

    BOL::BasicTimer runTimer;
    PotentialMinimizer& potentialMinimizer;
    double thresholdDistanceSquaredFraction;
    double extremumSeparationSquaredThreshold;
    PotentialMinimum dsbMinimum;
    PotentialMinimum globalMinimum;
    PotentialMinimum panicMinimum;
    bool dsbIsMetastable;
    int numberOfNudges;
    TunnelingCalculator& tunnelingCalculator;
    double quantumLifetimeThreshold;
    double thermalSurvivalThreshold;

    // SortMinima sets dsbMinimum to be the DSB minimum, panicMinimum to be the
    // minimum nearest the DSB minimum from the set of minima deeper than the
    // DSB minimum (or equal to dsbMinimum if there was no deeper minimum), and
    // globalMinimum to be the deepest minimum.
    void SortMinima();
    std::string XmlMinimum( PotentialMinimum const& potentialMinimum ) const;
  };

} /* namespace VevaciousPlusPlus */
#endif /* VEVACIOUSPLUSPLUS_HPP_ */
