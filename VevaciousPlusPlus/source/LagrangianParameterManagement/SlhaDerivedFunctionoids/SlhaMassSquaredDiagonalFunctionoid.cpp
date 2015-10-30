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
                             SlhaSourcedParameterFunctionoid const& squareMass,
                          SlhaSourcedParameterFunctionoid const& linearMass ) :
    SlhaDerivedFunctionoid(),
    squareMass( squareMass ),
    squareMassIndex( squareMass.IndexInValuesVector() ),
    linearMass( linearMass ),
    linearMassIndex( linearMass.IndexInValuesVector() )
  {
    // This constructor is just an initialization list.
  }

  SlhaMassSquaredDiagonalFunctionoid::~SlhaMassSquaredDiagonalFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
