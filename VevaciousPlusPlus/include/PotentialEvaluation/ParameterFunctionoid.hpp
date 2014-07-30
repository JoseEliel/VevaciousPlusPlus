/*
 * ParameterFunctionoid.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PARAMETERFUNCTIONOID_HPP_
#define PARAMETERFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class ParameterFunctionoid : public BOL::BasicObserved
  {
  public:
    ParameterFunctionoid( std::string const& creationString,
                          std::string const& pythonParameterName,
                          double const currentValue = 0.0 );
    virtual ~ParameterFunctionoid();


    // This returns the value of the functionoid based on the last update.
    double operator()() const{ return currentValue; }

    // This should return the value of the functionoid for the given logarithm
    // of the scale.
    virtual double operator()( double const logarithmOfScale ) const = 0;

    // This should update currentValue based on logarithmOfScale and any other
    // functionoids that this functionoid is dependent on.
    virtual void
    UpdateForNewLogarithmOfScale( double const logarithmOfScale ) = 0;

    // This is mainly for debugging.
    virtual std::string AsString() = 0;

    // This is mainly for debugging.
    std::string const& CreationString(){ return creationString; }

    // This is for creating a Python version of the potential.
    std::string const& PythonParameterName(){ return pythonParameterName; }

    // This is for creating a Python version of the potential.
    virtual std::string PythonParameterEvaluation() const = 0;


  protected:
    double currentValue;
    std::string const creationString;
    std::string const pythonParameterName;
  };

} /* namespace VevaciousPlusPlus */
#endif /* PARAMETERFUNCTIONOID_HPP_ */
