/*
 * BounceWithSplines.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEWITHSPLINES_HPP_
#define BOUNCEWITHSPLINES_HPP_

#include "../CommonIncludes.hpp"
#include "../PotentialEvaluation.hpp"
#include "../PotentialMinimization.hpp"
#include "TunnelingCalculator.hpp"

namespace VevaciousPlusPlus
{
  class BounceWithSplines : public TunnelingCalculator
  {
  public:
    BounceWithSplines( TunnelingStrategy const tunnelingStrategy,
                       double const survivalProbabilityThreshold );
    virtual
    ~BounceWithSplines();
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEWITHSPLINES_HPP_ */
