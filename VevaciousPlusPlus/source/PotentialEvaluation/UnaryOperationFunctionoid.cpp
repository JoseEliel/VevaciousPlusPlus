/*
 * UnaryOperationFunctionoid.cpp
 *
 *  Created on: Mar 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  UnaryOperationFunctionoid::UnaryOperationFunctionoid(
                                            double (*unaryOperation)( double ),
                             ParameterFunctionoid* const functionoidPointer ) :
    ParameterFunctionoid(),
    unaryOperation( unaryOperation ),
    functionoidPointer( functionoidPointer )
  {
    // This constructor is just an initialization list.
  }

  UnaryOperationFunctionoid::~UnaryOperationFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
