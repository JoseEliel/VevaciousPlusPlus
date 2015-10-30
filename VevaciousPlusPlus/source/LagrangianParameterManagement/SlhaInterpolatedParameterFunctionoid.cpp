/*
 * SlhaInterpolatedParameterFunctionoid.cpp
 *
 *  Created on: Oct 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaInterpolatedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  SlhaInterpolatedParameterFunctionoid::SlhaInterpolatedParameterFunctionoid(
                                              size_t const indexInValuesVector,
                              LHPC::SlhaSimplisticInterpreter const& lhaParser,
                                           std::string const& parameterName ) :
    SlhaSourcedParameterFunctionoid( indexInValuesVector ),
    parameterName( parameterName ),
    lhaParser( lhaParser )
  {
    // This constructor is just an initialization list.
  }

  SlhaInterpolatedParameterFunctionoid::SlhaInterpolatedParameterFunctionoid(
                     SlhaInterpolatedParameterFunctionoid const& copySource ) :
    SlhaSourcedParameterFunctionoid( copySource.indexInValuesVector ),
    parameterName( copySource.parameterName ),
    lhaParser( copySource.lhaParser )
  {
    // This constructor is just an initialization list.
  }

  SlhaInterpolatedParameterFunctionoid::~SlhaInterpolatedParameterFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
