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
            SlhaInterpolatedParameterFunctionoid const& firstChoiceFunctionoid,
        SlhaInterpolatedParameterFunctionoid const& secondChoiceFunctionoid ) :
    SlhaDerivedFunctionoid(),
    firstChoiceFunctionoid( firstChoiceFunctionoid ),
    firstChoiceIndex( firstChoiceFunctionoid.IndexInValuesVector() ),
    secondChoiceFunctionoid( secondChoiceFunctionoid ),
    secondChoiceIndex( secondChoiceFunctionoid.IndexInValuesVector() )
  {
    // This constructor is just an initialization list.
  }

  SlhaTwoSourceFunctionoid::~SlhaTwoSourceFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
