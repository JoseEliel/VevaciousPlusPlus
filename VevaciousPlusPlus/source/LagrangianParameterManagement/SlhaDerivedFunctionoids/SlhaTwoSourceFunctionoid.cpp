/*
 * SlhaTwoSourceFunctionoid.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaDerivedFunctionoids/SlhaTwoSourceFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  SlhaTwoSourceFunctionoid::SlhaTwoSourceFunctionoid(
                                              size_t const indexInValuesVector,
                                                 size_t const firstChoiceIndex,
                                            size_t const secondChoiceIndex  ) :
    SlhaSourcedParameterFunctionoid( indexInValuesVector ),
    firstChoiceIndex( firstChoiceIndex ),
    secondChoiceIndex( secondChoiceIndex )
  {
    // This constructor is just an initialization list.
  }

  SlhaTwoSourceFunctionoid::~SlhaTwoSourceFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
