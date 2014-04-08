/*
 * TunnelingCalculator.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef TUNNELINGCALCULATOR_HPP_
#define TUNNELINGCALCULATOR_HPP_

#include "../StandardIncludes.hpp"
#include "../PotentialEvaluation/PotentialFunction.hpp"

namespace VevaciousPlusPlus
{
  class TunnelingCalculator
  {
  public:
    TunnelingCalculator( PotentialFunction& potentialFunction );
    virtual
    ~TunnelingCalculator();


  protected:
    PotentialFunction& potentialFunction;
  };

} /* namespace VevaciousPlusPlus */
#endif /* TUNNELINGCALCULATOR_HPP_ */
