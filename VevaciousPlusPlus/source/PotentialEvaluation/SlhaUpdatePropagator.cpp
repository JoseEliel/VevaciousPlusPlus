/*
 * SlhaUpdatePropagator.cpp
 *
 *  Created on: Apr 17, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  SlhaUpdatePropagator::SlhaUpdatePropagator(
                                   SlhaUpdatePropagator& previousPropagator ) :
    BOL::PushedToObserver< SlhaManager >(),
    BOL::PushingObserved< SlhaManager >(),
    slhaManager( previousPropagator.slhaManager )
  {
    previousPropagator.registerObserver( this );
  }

  SlhaUpdatePropagator::SlhaUpdatePropagator( SlhaManager& slhaManager ) :
    BOL::PushedToObserver< SlhaManager >(),
    BOL::PushingObserved< SlhaManager >(),
    slhaManager( slhaManager )
  {
    slhaManager.registerObserver( this );
  }

  SlhaUpdatePropagator::~SlhaUpdatePropagator()
  {
    // This does nothing beyond the destructors of the base classes.
  }

} /* namespace VevaciousPlusPlus */
