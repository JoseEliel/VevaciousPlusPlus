/*
 * ParameterFunctionoid.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef RUNNINGPARAMETER_HPP_
#define RUNNINGPARAMETER_HPP_

#include "../StandardIncludes.hpp"

namespace VevaciousPlusPlus
{

  class ParameterFunctionoid
  {
  public:
    ParameterFunctionoid( double const currentValue = 0.0 );
    virtual
    ~ParameterFunctionoid();


    // This returns the value of the functionoid based on the last update.
    double operator()() const{ return currentValue; }

    // This should update currentValue based on logarithmOfScale and any other
    // functionoids that this functionoid is dependent on.
    virtual void
    UpdateForNewLogarithmOfScale( double const logarithmOfScale ) = 0;

    // This is mainly for debugging.
    virtual std::string AsString() = 0;


  protected:
    double currentValue;
  };

} /* namespace VevaciousPlusPlus */
#endif /* RUNNINGPARAMETER_HPP_ */
