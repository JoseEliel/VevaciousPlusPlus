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
                        SlhaSourcedParameterFunctionoid const& directTrilinear,
                    SlhaSourcedParameterFunctionoid const& trilinearOverYukawa,
                   SlhaSourcedParameterFunctionoid const& appropriateYukawa ) :
    SlhaDerivedFunctionoid(),
    directTrilinear( directTrilinear ),
    directTrilinearIndex( directTrilinear.IndexInValuesVector() ),
    trilinearOverYukawa( trilinearOverYukawa ),
    trilinearOverYukawaIndex( trilinearOverYukawa.IndexInValuesVector() ),
    appropriateYukawa( appropriateYukawa ),
    appropriateYukawaIndex( appropriateYukawa.IndexInValuesVector() )
  {
    // This constructor is just an initialization list.
  }

  SlhaTrilinearDiagonalFunctionoid::~SlhaTrilinearDiagonalFunctionoid()
  {
    // This does nothing.s
  }

} /* namespace VevaciousPlusPlus */
