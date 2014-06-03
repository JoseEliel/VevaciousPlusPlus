/*
 * OdeintBubbleObserver.cpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{

  OdeintBubbleObserver::OdeintBubbleObserver(
      std::vector< BubbleRadialValueDescription >& bubbleDescription ) :
    bubbleDescription( bubbleDescription ),
    definitelyUndershot( false ),
    definitelyOvershot( false ),
    overshootIndex( 0 ),
    needsOrdering( false )
  {
    // This constructor is just an initialization list.
  }

  OdeintBubbleObserver::~OdeintBubbleObserver()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
