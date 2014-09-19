/*
 * OdeintBubbleObserver.cpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/OdeintBubbleObserver.hpp"

namespace VevaciousPlusPlus
{

  OdeintBubbleObserver::OdeintBubbleObserver(
             std::vector< BubbleRadialValueDescription >& bubbleDescription ) :
    bubbleDescription( bubbleDescription )
  {
    // This constructor is just an initialization list.
  }

  OdeintBubbleObserver::~OdeintBubbleObserver()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
