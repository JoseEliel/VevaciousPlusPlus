/*
 * UnaryOperationFunctionoid.hpp
 *
 *  Created on: Mar 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef UNARYOPERATIONFUNCTIONOID_HPP_
#define UNARYOPERATIONFUNCTIONOID_HPP_

#include "../StandardIncludes.hpp"
#include "ParameterFunctionoid.hpp"


namespace VevaciousPlusPlus
{

  class UnaryOperationFunctionoid : public ParameterFunctionoid
  {
  public:
    UnaryOperationFunctionoid( double (*unaryOperation)( double ),
                              ParameterFunctionoid* const functionoidPointer );
    virtual
    ~UnaryOperationFunctionoid();


    // This re-calculates the coefficients of the polynomial of the logarithm
    // of the scale used in evaluating the functionoid.
    virtual void UpdateForNewLogarithmOfScale( double const logarithmOfScale );

    // This is mainly for debugging.
    virtual std::string AsString();


  protected:
    double (*unaryOperation)( double );
    ParameterFunctionoid* const functionoidPointer;
  };




  // This updates currentValue based on logarithmOfScale.
  inline void UnaryOperationFunctionoid::UpdateForNewLogarithmOfScale(
                                                double const logarithmOfScale )
  {
    currentValue = (*unaryOperation)( (*functionoidPointer)() );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "[" << this->AsString()
    << "].UpdateForNewLogarithmOfScale( " << logarithmOfScale
    << " ) called. currentValue = " << currentValue;
    std::cout << std::endl;*/
  }

  // This is mainly for debugging.
  inline std::string UnaryOperationFunctionoid::AsString()
  {
    std::stringstream returnStream;
    returnStream << "[UNARYOPERATION " << this << ": "
    << functionoidPointer->AsString() << "]";
    return std::string( returnStream.str() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* UNARYOPERATIONFUNCTIONOID_HPP_ */
