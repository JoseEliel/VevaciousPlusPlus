/*
 * VevaciousPlusPlus.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef VEVACIOUSPLUSPLUS_HPP_
#define VEVACIOUSPLUSPLUS_HPP_

#include "CommonIncludes.hpp"
#include "PotentialEvaluation/PotentialFunctions/PotentialFunction.hpp"
#include "PotentialEvaluation/PotentialFunctions/PotentialFromPolynomialAndMasses.hpp"
#include "PotentialMinimization/PotentialMinimizer.hpp"
#include "TunnelingCalculation/TunnelingCalculator.hpp"

namespace VevaciousPlusPlus
{
  class VevaciousPlusPlus
  {
  public:
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
