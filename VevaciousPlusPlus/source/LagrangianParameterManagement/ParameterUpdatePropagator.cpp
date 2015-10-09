/*
 * SlhaUpdatePropagator.cpp
 *
 *  Created on: Apr 17, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/LagrangianParameterManagement/ParameterUpdatePropagator.hpp"

namespace VevaciousPlusPlus
{

  ParameterUpdatePropagator::ParameterUpdatePropagator(
                                   ParameterUpdatePropagator& previousPropagator ) :
    BOL::PushedToObserver< SlhaManager >(),
    BOL::PushingObserved< SlhaManager >(),
    slhaManager( previousPropagator.slhaManager )
  {
    previousPropagator.registerObserver( this );
  }

  ParameterUpdatePropagator::ParameterUpdatePropagator( SlhaManager& slhaManager ) :
    BOL::PushedToObserver< SlhaManager >(),
    BOL::PushingObserved< SlhaManager >(),
    slhaManager( slhaManager )
  {
    slhaManager.registerObserver( this );
  }

  ParameterUpdatePropagator::~ParameterUpdatePropagator()
  {
    // This does nothing beyond the destructors of the base classes.
  }

} /* namespace VevaciousPlusPlus */
