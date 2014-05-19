/*
 * OdeintBubbleObserver.cpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{

  OdeintBubbleObserver::OdeintBubbleObserver() :
    radialValues(),
    auxiliaryValues(),
    auxiliaryDerivatives(),
    definitelyUndershot( false ),
    definitelyOvershot( false )
  {
    // This constructor is just an initialization list.
  }

  OdeintBubbleObserver::~OdeintBubbleObserver()
  {
    // This does nothing.
  }


  // This is in the form required for the Boost odeint package.
  void OdeintBubbleObserver::operator()(
                      std::vector< double > const& auxiliaryAndFirstDerivative,
                                                double const radialValue )
  {
    size_t insertionOffset( radialValues.size() );
    while( ( insertionOffset > 0 )
           &&
           ( radialValues[ insertionOffset - 1 ] > radialValue ) )
    {
      --insertionOffset;
    }
    radialValues.insert( ( radialValues.begin() + insertionOffset ),
                         radialValue );
    auxiliaryValues.insert( ( auxiliaryValues.begin() + insertionOffset ),
                            auxiliaryAndFirstDerivative[ 0 ] );
    auxiliaryDerivatives.insert( ( auxiliaryDerivatives.begin()
                                   + insertionOffset ),
                                 auxiliaryAndFirstDerivative[ 1 ] );
    // We update definitelyUndershot and definitelyOvershot if they're not yet
    // true (because whether or not there was an undershoot or overshoot does
    // not only depend on the last timestep, but on whether any of the
    // timesteps produced a < 0.0 or da/dr > 0.0).
    if( !definitelyUndershot )
    {
      definitelyUndershot = ( auxiliaryAndFirstDerivative[ 1 ] > 0.0 );
    }
    if( !definitelyOvershot )
    {
      definitelyOvershot = ( auxiliaryAndFirstDerivative[ 0 ] < 0.0 );
    }
  }

} /* namespace VevaciousPlusPlus */
