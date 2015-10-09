/*
 * ParameterUpdatePropagator.hpp
 *
 *  Created on: Apr 17, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PARAMETERUPDATEPROPAGATOR_HPP_
#define PARAMETERUPDATEPROPAGATOR_HPP_

#include "../LagrangianParameterManagement/LagrangianParameterManager.hpp"
#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class ParameterUpdatePropagator :
                    public BOL::PushedToObserver< LagrangianParameterManager >,
                      public BOL::PushingObserved< LagrangianParameterManager >
  {
  public:
    ParameterUpdatePropagator( ParameterUpdatePropagator& previousPropagator );
    ParameterUpdatePropagator(
                      LagrangianParameterManager& lagrangianParameterManager );
    virtual ~ParameterUpdatePropagator();


    LagrangianParameterManager const& GetLagrangianParameterManager() const
    { return lagrangianParameterManager; }

    // This should perform all relevant updates for the new SLHA data except
    // for propagating the push to the set of dependent SlhaUpdatePropagators.
    virtual void UpdateSelfForNewParameterPoint(
            LagrangianParameterManager const& lagrangianParameterManager ) = 0;

    // This calls the methods first to update this instance and then push
    // lagrangianParameterManager to the ParameterUpdatePropagator objects
    // observing this one.
    virtual void respondToPush(
                LagrangianParameterManager const& lagrangianParameterManager );

    virtual void respondToObservedSignal()
    { respondToPush( lagrangianParameterManager ); }


  protected:
    LagrangianParameterManager& lagrangianParameterManager;
  };





  // This calls the methods first to update this instance and then push
  // lagrangianParameterManager to the ParameterUpdatePropagator objects
  // observing this one.
  inline void ParameterUpdatePropagator::respondToPush(
                 LagrangianParameterManager const& lagrangianParameterManager )
  {
    UpdateSelfForNewParameterPoint( lagrangianParameterManager );
    updateObservers( lagrangianParameterManager );
  }

} /* namespace VevaciousPlusPlus */
#endif /* PARAMETERUPDATEPROPAGATOR_HPP_ */
