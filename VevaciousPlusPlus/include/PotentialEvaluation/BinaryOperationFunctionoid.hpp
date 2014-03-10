/*
 * BinaryOperationFunctionoid.hpp
 *
 *  Created on: Mar 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PLUSMINUSFUNCTIONOID_HPP_
#define PLUSMINUSFUNCTIONOID_HPP_

#include "../StandardIncludes.hpp"
#include "ParameterFunctionoid.hpp"


namespace VevaciousPlusPlus
{

  class BinaryOperationFunctionoid : public ParameterFunctionoid
  {
  public:
    static double PlusFunction( double firstValue,
                                double secondValue )
    { return ( firstValue + secondValue ); }
    static double MinusFunction( double firstValue,
                                 double secondValue )
    { return ( firstValue - secondValue ); }
    static double MultiplyFunction( double firstValue,
                                    double secondValue )
    { return ( firstValue * secondValue ); }
    static double DivideFunction( double firstValue,
                                  double secondValue )
    { return ( firstValue / secondValue ); }

    BinaryOperationFunctionoid( double (*binaryOperation)( double,
                                                           double ),
                                ParameterFunctionoid* const firstFunctionoid,
                               ParameterFunctionoid* const secondFunctionoid );
    virtual
    ~BinaryOperationFunctionoid();


    // This re-calculates the coefficients of the polynomial of the logarithm
    // of the scale used in evaluating the functionoid.
    void UpdateForNewLogarithmOfScale( double const logarithmOfScale );

    // This is mainly for debugging.
    virtual std::string AsString();


  protected:
    double (*binaryOperation)( double,
                               double );
    ParameterFunctionoid* firstFunctionoid;
    ParameterFunctionoid* secondFunctionoid;
  };




  // This updates currentValue based on logarithmOfScale.
  inline void BinaryOperationFunctionoid::UpdateForNewLogarithmOfScale(
                                                double const logarithmOfScale )
  {
    currentValue = (*binaryOperation)( (*firstFunctionoid)(),
                                       (*secondFunctionoid)() );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "[" << this->AsString()
    << "].UpdateForNewLogarithmOfScale( " << logarithmOfScale
    << " ) called. currentValue = " << currentValue;
    std::cout << std::endl;*/
  }

  // This is mainly for debugging.
  inline std::string BinaryOperationFunctionoid::AsString()
  {
    std::stringstream returnStream;
    returnStream
    << "[BINARYOPERATION " << this << ": " << firstFunctionoid->AsString()
    << ", " << secondFunctionoid->AsString() << "]";
    return std::string( returnStream.str() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* PLUSMINUSFUNCTIONOID_HPP_ */
