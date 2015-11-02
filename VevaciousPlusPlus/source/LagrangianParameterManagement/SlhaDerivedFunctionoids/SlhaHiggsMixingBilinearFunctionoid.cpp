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
            SlhaSourcedParameterFunctionoid const& treePseudoscalarMassSquared,
                             SlhaSourcedParameterFunctionoid const& tanBeta ) :
    SlhaSourcedParameterFunctionoid( indexInValuesVector ),
    treePseudoscalarMassSquared( treePseudoscalarMassSquared ),
    treePseudoscalarMassSquaredIndex(
                           treePseudoscalarMassSquared.IndexInValuesVector() ),
    tanBeta( tanBeta ),
    tanBetaIndex( tanBeta.IndexInValuesVector() )
  {
    // This constructor is just an initialization list.
  }

  SlhaHiggsMixingBilinearFunctionoid::~SlhaHiggsMixingBilinearFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
