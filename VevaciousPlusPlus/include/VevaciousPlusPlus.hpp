/*
 * VevaciousPlusPlus.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef VEVACIOUSPLUSPLUS_HPP_
#define VEVACIOUSPLUSPLUS_HPP_

#include "CommonIncludes.hpp"
#include "PotentialEvaluation.hpp"
#include "PotentialMinimization.hpp"
#include "TunnelingCalculation.hpp"
#include <ctime>

namespace VevaciousPlusPlus
{
  class VevaciousPlusPlus
  {
  public:
    static std::string const versionString;
    static std::string const citationString;

    VevaciousPlusPlus( std::string const& initializationFileName );

    VevaciousPlusPlus( BOL::ArgumentParser& argumentParser,
                       SlhaManager& slhaManager,
                       PotentialMinimizer& potentialMinimizer,
                       TunnelingCalculator& tunnelingCalculator );

    virtual ~VevaciousPlusPlus();

    void RunPoint( std::string const& parameterFilename );

    void WriteXmlResults( std::string const& xmlFilename );

    void WriteSlhaResults( std::string const& slhaFilename,
                           bool const writeWarnings = true );


  protected:
    PotentialFromPolynomialAndMasses* potentialFunction;
    PotentialFromPolynomialAndMasses* potentialFunctionDeleter;
    RunningParameterManager runningParameterManager;
    PotentialMinimizer* potentialMinimizer;
    PotentialMinimizer* potentialMinimizerDeleter;
    TunnelingCalculator* tunnelingCalculator;
    TunnelingCalculator* tunnelingCalculatorDeleter;
    time_t currentTime;
  };

} /* namespace VevaciousPlusPlus */
#endif /* VEVACIOUSPLUSPLUS_HPP_ */
