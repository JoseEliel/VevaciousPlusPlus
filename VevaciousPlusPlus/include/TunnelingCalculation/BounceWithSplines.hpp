/*
 * BounceWithSplines.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEWITHSPLINES_HPP_
#define BOUNCEWITHSPLINES_HPP_

#include "../StandardIncludes.hpp"
#include "../PotentialEvaluation/PotentialEvaluation.hpp"
#include "TunnelingCalculator.hpp"

namespace VevaciousPlusPlus
{
  class BounceWithSplines : public TunnelingCalculator
  {
  public:
    BounceWithSplines( PotentialFunction& potentialFunction );
    virtual
    ~BounceWithSplines();
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEWITHSPLINES_HPP_ */
