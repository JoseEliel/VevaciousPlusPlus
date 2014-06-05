/*
 * OdeintBubbleObserver.hpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef ODEINTBUBBLEOBSERVER_HPP_
#define ODEINTBUBBLEOBSERVER_HPP_

#include "../CommonIncludes.hpp"
#include "BubbleRadialValueDescription.hpp"

namespace VevaciousPlusPlus
{

  class OdeintBubbleObserver
  {
  public:
    OdeintBubbleObserver(
              std::vector< BubbleRadialValueDescription >& bubbleDescription );
    virtual
    ~OdeintBubbleObserver();


    // This pushes back the radial value, the auxiliary value, and its slope
    // into bubbleDescription.
    void operator()( std::vector< double > const& auxiliaryAndFirstDerivative,
                     double const radialValue )
    { bubbleDescription.push_back( BubbleRadialValueDescription( radialValue,
                                              auxiliaryAndFirstDerivative[ 0 ],
                                        auxiliaryAndFirstDerivative[ 1 ] ) );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "OdeintBubbleObserver::operator() finished: bubbleDescription.size() = "
    << bubbleDescription.size() << ", bubbleDescription.back() = { r = "
    << bubbleDescription.back().radialValue << ", p = "
    << bubbleDescription.back().auxiliaryValue << ", dp/dr = "
    << bubbleDescription.back().auxiliarySlope << " }";
    std::cout << std::endl;/**/ }


  protected:
    std::vector< BubbleRadialValueDescription >& bubbleDescription;
  };

} /* namespace VevaciousPlusPlus */
#endif /* ODEINTBUBBLEOBSERVER_HPP_ */
