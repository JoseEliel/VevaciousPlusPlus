/*
 * GradientMinimizer.cpp
 *
 *  Created on: Jun 30, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/GradientMinimizer.hpp"

namespace VevaciousPlusPlus
{

  GradientMinimizer::GradientMinimizer(
                                 PotentialFunction const& potentialFunction ) :
    potentialFunction( potentialFunction )
  {
    // This constructor is just an initialization list.
  }

  GradientMinimizer::~GradientMinimizer()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
