/*
 * SlhaMassSquaredDiagonalFunctionoid.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaDerivedFunctionoids/SlhaMassSquaredDiagonalFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  SlhaMassSquaredDiagonalFunctionoid::SlhaMassSquaredDiagonalFunctionoid(
                                              size_t const indexInValuesVector,
                                                  size_t const squareMassIndex,
                                               size_t const linearMassIndex ) :
    LhaSourcedParameterFunctionoid( indexInValuesVector ),
    squareMassIndex( squareMassIndex ),
    linearMassIndex( linearMassIndex )
  {
    // This constructor is just an initialization list.
  }

  SlhaMassSquaredDiagonalFunctionoid::~SlhaMassSquaredDiagonalFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
