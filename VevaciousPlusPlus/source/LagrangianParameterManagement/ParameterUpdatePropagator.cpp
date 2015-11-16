/*
 * SlhaUpdatePropagator.cpp
 *
 *  Created on: Apr 17, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/ParameterUpdatePropagator.hpp"

namespace VevaciousPlusPlus
{

  ParameterUpdatePropagator::ParameterUpdatePropagator(
                              ParameterUpdatePropagator& previousPropagator ) :
    BOL::PushedToObserver< LagrangianParameterManager >(),
    BOL::PushingObserved< LagrangianParameterManager >(),
    lagrangianParameterManager( previousPropagator.lagrangianParameterManager )
  {
    previousPropagator.registerObserver( this );
  }

  ParameterUpdatePropagator::ParameterUpdatePropagator(
                     LagrangianParameterManager& lagrangianParameterManager ) :
    BOL::PushedToObserver< LagrangianParameterManager >(),
    BOL::PushingObserved< LagrangianParameterManager >(),
    lagrangianParameterManager( lagrangianParameterManager )
  {
    lagrangianParameterManager.registerObserver( this );
  }

  ParameterUpdatePropagator::~ParameterUpdatePropagator()
  {
    // This does nothing beyond the destructors of the base classes.
  }

} /* namespace VevaciousPlusPlus */
