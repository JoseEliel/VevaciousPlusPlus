/*
 * HomotopyContinuationReadyPotential.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONREADYPOTENTIAL_HPP_
#define HOMOTOPYCONTINUATIONREADYPOTENTIAL_HPP_

#include "../StandardIncludes.hpp"
#include "PotentialFunction.hpp"

namespace VevaciousPlusPlus
{
  class HomotopyContinuationReadyPotential : public PotentialFunction
  {
  public:
    HomotopyContinuationReadyPotential();
    virtual
    ~HomotopyContinuationReadyPotential();
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONREADYPOTENTIAL_HPP_ */
