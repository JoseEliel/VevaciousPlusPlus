/*
 * HomotopyContinuationAndGradient.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONANDGRADIENT_HPP_
#define HOMOTOPYCONTINUATIONANDGRADIENT_HPP_

#include "../StandardIncludes.hpp"
#include "PotentialMinimizer.hpp"
#include "../PotentialEvaluation/HomotopyContinuationReadyPotential.hpp"

namespace VevaciousPlusPlus
{
  class HomotopyContinuationAndGradient : public PotentialMinimizer
  {
  public:
    HomotopyContinuationAndGradient(
                       HomotopyContinuationReadyPotential& potentialFunction );
    virtual
    ~HomotopyContinuationAndGradient();
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONANDGRADIENT_HPP_ */
