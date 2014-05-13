/*
 * BubbleRadiusFromAuxiliary.cpp
 *
 *  Created on: May 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  BubbleRadiusFromAuxiliary::BubbleRadiusFromAuxiliary(
                                         double const appropriateScaleInGev ) :
    appropriateScaleInInverseGev( 1.0 / appropriateScaleInGev )
  {
    // This constructor is just an initialization list.
  }

  BubbleRadiusFromAuxiliary::~BubbleRadiusFromAuxiliary()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
