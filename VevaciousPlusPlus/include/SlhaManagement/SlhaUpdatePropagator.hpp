/*
 * SlhaUpdatePropagator.hpp
 *
 *  Created on: Apr 17, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHAUPDATEPROPAGATOR_HPP_
#define SLHAUPDATEPROPAGATOR_HPP_

#include "CommonIncludes.hpp"
#include "SlhaManager.hpp"

namespace VevaciousPlusPlus
{

  class SlhaUpdatePropagator : public BOL::PushedToObserver< SlhaManager >,
                               public BOL::PushingObserved< SlhaManager >
  {
  public:
    SlhaUpdatePropagator( SlhaUpdatePropagator& previousPropagator );
    SlhaUpdatePropagator( SlhaManager& slhaManager );
    virtual ~SlhaUpdatePropagator();


    // This should perform all relevant updates for the new SLHA data except
    // for propagating the push to the set of dependent SlhaUpdatePropagators.
    virtual void UpdateSelfForNewSlha( SlhaManager const& slhaManager ) = 0;

    virtual void respondToPush( SlhaManager const& slhaManager )
    { UpdateSelfForNewSlha( slhaManager ); updateObservers( slhaManager ); }

    virtual void respondToObservedSignal(){ respondToPush( slhaManager ); }

  protected:
    SlhaManager& slhaManager;
  };

} /* namespace VevaciousPlusPlus */
#endif /* SLHAUPDATEPROPAGATOR_HPP_ */
