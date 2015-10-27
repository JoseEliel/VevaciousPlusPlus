/*
 * SlhaTanBetaFunctionoid.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaDerivedFunctionoids/SlhaTanBetaFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  SlhaTanBetaFunctionoid::SlhaTanBetaFunctionoid(
                       SlhaInterpolatedParameterFunctionoid const& hmixTanBeta,
                  SlhaInterpolatedParameterFunctionoid const& extparTanBeta ) :
    hmixTanBeta( hmixTanBeta ),
    hmixTanBetaIndex( hmixTanBeta.IndexInValuesVector() ),
    extparTanBeta( extparTanBeta ),
    extparTanBetaIndex( extparTanBeta.IndexInValuesVector() )
  {
    // This constructor is just an initialization list.
  }

  SlhaTanBetaFunctionoid::~SlhaTanBetaFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
