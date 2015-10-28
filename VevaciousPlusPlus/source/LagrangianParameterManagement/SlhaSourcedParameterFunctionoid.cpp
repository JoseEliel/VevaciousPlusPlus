/*
 * SlhaSourcedParameterFunctionoid.cpp
 *
 *  Created on: Oct 28, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaSourcedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  SlhaSourcedParameterFunctionoid::SlhaSourcedParameterFunctionoid() :
    indexInValuesVector( -1 )
  {
    // This constructor is just an initialization list.
  }

  SlhaSourcedParameterFunctionoid::SlhaSourcedParameterFunctionoid(
                                           size_t const indexInValuesVector ) :
    indexInValuesVector( indexInValuesVector )
  {
    // This constructor is just an initialization list.
  }

  SlhaSourcedParameterFunctionoid::~SlhaSourcedParameterFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
