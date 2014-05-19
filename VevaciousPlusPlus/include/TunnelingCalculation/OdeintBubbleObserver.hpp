/*
 * OdeintBubbleObserver.hpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef ODEINTBUBBLEOBSERVER_HPP_
#define ODEINTBUBBLEOBSERVER_HPP_

#include "../CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class OdeintBubbleObserver
  {
  public:
    OdeintBubbleObserver();
    virtual
    ~OdeintBubbleObserver();


    // This is in the form required for the Boost odeint package.
    void operator()( std::vector< double > const& auxiliaryAndFirstDerivative,
                     double const radialValue );

    std::vector< double > const& RadialValues() const{ return radialValues; }
    std::vector< double > const&
    AuxiliaryValues() const{ return auxiliaryValues; }
    std::vector< double > const&
    AuxiliaryDerivatives() const{ return auxiliaryDerivatives; }

    // This returns definitelyUndershot which is set true if the slope of the
    // auxiliary variable goes positive for any radial value.
    bool DefinitelyUndershot() const{ return definitelyUndershot; }

    // This returns definitelyOvershot which is set true if the auxiliary
    // variable goes negative for any radial value.
    bool DefinitelyOvershot() const{ return definitelyOvershot; }

    // This resets everything for a new integration.
    void ResetValues();


  protected:
    std::vector< double > radialValues;
    std::vector< double > auxiliaryValues;
    std::vector< double > auxiliaryDerivatives;
    bool definitelyUndershot;
    bool definitelyOvershot;
  };

  // This resets everything for a new integration.
  inline void OdeintBubbleObserver::ResetValues()
  {
    radialValues.clear();
    auxiliaryValues.clear();
    auxiliaryDerivatives.clear();
    definitelyUndershot = false;
    definitelyOvershot = false;
  }

} /* namespace VevaciousPlusPlus */
#endif /* ODEINTBUBBLEOBSERVER_HPP_ */
