/*
 * CosmoTransitionsRunner.hpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef COSMOTRANSITIONSRUNNER_HPP_
#define COSMOTRANSITIONSRUNNER_HPP_

#include "../CommonIncludes.hpp"
#include "../PotentialEvaluation.hpp"
#include "../PotentialMinimization.hpp"
#include "BounceWithSplines.hpp"


namespace VevaciousPlusPlus
{

  class CosmoTransitionsRunner : public BounceWithSplines
  {
  public:
    CosmoTransitionsRunner( IWritesPythonPotential& potentialFunction,
                            TunnelingStrategy const tunnelingStrategy,
                            double const survivalProbabilityThreshold,
                            std::string const& pathToCosmotransitions );
    virtual
    ~CosmoTransitionsRunner();


    // This creates an entire Python program to run CosmoTransitions to
    // calculate the bounce action(s) for the tunneling strategy.
    virtual void CalculateTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum );


  protected:
    IWritesPythonPotential& potentialFunction;
    std::string const pathToCosmotransitions;
  };

} /* namespace VevaciousPlusPlus */
#endif /* COSMOTRANSITIONSRUNNER_HPP_ */
