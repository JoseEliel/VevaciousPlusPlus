/*
 * SlhaTrilinearDiagonalFunctionoid.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaDerivedFunctionoids/SlhaTrilinearDiagonalFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  SlhaTrilinearDiagonalFunctionoid::SlhaTrilinearDiagonalFunctionoid(
                                              size_t const indexInValuesVector,
                                             size_t const directTrilinearIndex,
                                         size_t const trilinearOverYukawaIndex,
                                        size_t const appropriateYukawaIndex ) :
    LhaSourcedParameterFunctionoid( indexInValuesVector ),
    directTrilinearIndex( directTrilinearIndex ),
    trilinearOverYukawaIndex( trilinearOverYukawaIndex ),
    appropriateYukawaIndex( appropriateYukawaIndex )
  {
    // This constructor is just an initialization list.
  }

  SlhaTrilinearDiagonalFunctionoid::~SlhaTrilinearDiagonalFunctionoid()
  {
    // This does nothing.s
  }

} /* namespace VevaciousPlusPlus */
