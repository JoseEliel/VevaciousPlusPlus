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
                 LHPC::SLHA::SparseManyIndexedBlock< double > const& slhaBlock,
                                             std::string const& indexString ) :
    SlhaSourcedParameterFunctionoid( indexInValuesVector ),
    slhaBlock( slhaBlock ),
    indexVector( BOL::StringParser::stringToIntVector( indexString ) )
  {
    // This constructor is just an initialization list.
  }

  SlhaInterpolatedParameterFunctionoid::SlhaInterpolatedParameterFunctionoid(
                     SlhaInterpolatedParameterFunctionoid const& copySource ) :
    SlhaSourcedParameterFunctionoid( copySource.indexInValuesVector ),
    slhaBlock( copySource.slhaBlock ),
    indexVector( copySource.indexVector )
  {
    // This constructor is just an initialization list.
  }

  SlhaInterpolatedParameterFunctionoid::~SlhaInterpolatedParameterFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
