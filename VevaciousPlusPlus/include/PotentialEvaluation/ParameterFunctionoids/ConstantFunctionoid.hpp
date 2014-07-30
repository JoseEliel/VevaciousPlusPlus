/*
 * ConstantFunctionoid.hpp
 *
 *  Created on: Mar 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef CONSTANTFUNCTIONOID_HPP_
#define CONSTANTFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "../ParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class ConstantFunctionoid : public ParameterFunctionoid
  {
  public:
    ConstantFunctionoid( double const constantValue,
                         std::string const& creationString,
                         std::string const& pythonParameterName );
    virtual ~ConstantFunctionoid();


    // This returns the value of the functionoid for the given logarithm of the
    // scale, which is just the value it always had.
    virtual double operator()( double const logarithmOfScale ) const
    { return currentValue; }

    // This does nothing, actually.
    virtual void
    UpdateForNewLogarithmOfScale( double const logarithmOfScale ){}

    // This is mainly for debugging.
    virtual std::string AsString();

    // This is for creating a Python version of the potential.
    virtual std::string PythonParameterEvaluation() const;
  };




  // This is mainly for debugging.
  inline std::string ConstantFunctionoid::AsString()
  {
    std::stringstream returnStream;
    returnStream << "[CONSTANT " << currentValue << "]";
    return returnStream.str();
  }

  // This is for creating a Python version of the potential.
  inline std::string ConstantFunctionoid::PythonParameterEvaluation() const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << std::setprecision( 12 ) << pythonParameterName << " = " << currentValue;
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* CONSTANTFUNCTIONOID_HPP_ */
