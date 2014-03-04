/*
 * BinaryOperationFunctionoid.cpp
 *
 *  Created on: Mar 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  BinaryOperationFunctionoid::BinaryOperationFunctionoid(
                                       double (*binaryOperation)( double const,
                                                                double const ),
                                  ParameterFunctionoid* const firstFunctionoid,
                              ParameterFunctionoid* const secondFunctionoid ) :
    ParameterFunctionoid(),
    binaryOperation( binaryOperation ),
    firstFunctionoid( firstFunctionoid ),
    secondFunctionoid( secondFunctionoid )
  {
    // This constructor is just an initialization list.
  }

  BinaryOperationFunctionoid::~BinaryOperationFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
