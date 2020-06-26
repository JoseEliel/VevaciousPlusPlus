/*
 * OdeintBubbleObserver.hpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef ODEINTBUBBLEOBSERVER_HPP_
#define ODEINTBUBBLEOBSERVER_HPP_

#include <vector>
#include "BubbleRadialValueDescription.hpp"

namespace VevaciousPlusPlus
{

  class OdeintBubbleObserver
  {
  public:
    OdeintBubbleObserver(
             std::vector< BubbleRadialValueDescription >& bubbleDescription ) :
      bubbleDescription( bubbleDescription ) {}

    virtual ~OdeintBubbleObserver() {}


    // This pushes back the radial value, the auxiliary value, and its slope
    // into bubbleDescription.
    void operator()( std::vector< double > const& auxiliaryAndFirstDerivative,
                     double const radialValue )
    { bubbleDescription.push_back( BubbleRadialValueDescription( radialValue,
                                              auxiliaryAndFirstDerivative[ 0 ],
                                        auxiliaryAndFirstDerivative[ 1 ] ) ); }


  protected:
    std::vector< BubbleRadialValueDescription >& bubbleDescription;
  };

} /* namespace VevaciousPlusPlus */
#endif /* ODEINTBUBBLEOBSERVER_HPP_ */
