/*
 * SlhaHiggsMixingBilinearFunctionoid.cpp
 *
 *  Created on: Oct 28, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaDerivedFunctionoids/SlhaHiggsMixingBilinearFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  SlhaHiggsMixingBilinearFunctionoid::SlhaHiggsMixingBilinearFunctionoid(
                                              size_t const indexInValuesVector,
                                 size_t const treePseudoscalarMassSquaredIndex,
                                                  size_t const tanBetaIndex ) :
    SlhaSourcedParameterFunctionoid( indexInValuesVector ),
    treePseudoscalarMassSquaredIndex( treePseudoscalarMassSquaredIndex ),
    tanBetaIndex( tanBetaIndex )
  {
    // This constructor is just an initialization list.
  }

  SlhaHiggsMixingBilinearFunctionoid::~SlhaHiggsMixingBilinearFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
