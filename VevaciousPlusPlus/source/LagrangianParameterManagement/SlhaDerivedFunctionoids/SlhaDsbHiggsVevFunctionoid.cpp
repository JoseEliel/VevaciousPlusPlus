/*
 * SlhaDsbVdFunctionoid.cpp
 *
 *  Created on: Oct 26, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaDerivedFunctionoids/SlhaDsbHiggsVevFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  SlhaDsbHiggsVevFunctionoid::SlhaDsbHiggsVevFunctionoid(
            SlhaInterpolatedParameterFunctionoid const& sarahHiggsVevComponent,
       SlhaInterpolatedParameterFunctionoid const& slhaHiggsVevEuclideanLength,
                                         SlhaTanBetaFunctionoid const& tanBeta,
                                                       bool const sinNotCos ) :
    SlhaDerivedFunctionoid(),
    sarahHiggsVevComponent( sarahHiggsVevComponent ),
    sarahHmixIndex( sarahHiggsVevComponent.IndexInValuesVector() ),
    slhaHiggsVevEuclideanLength( slhaHiggsVevEuclideanLength ),
    slhaHmixIndex( slhaHiggsVevEuclideanLength.IndexInValuesVector() ),
    tanBeta( tanBeta ),
    tanBetaIndex( tanBeta.IndexInValuesVector() ),
    cosOrSin( &cos ),
    cosOrSinPythonString( "math.cos" )
  {
    // The assumption is that the functionoid is for vd, hence uses cos(beta).
    if( sinNotCos )
    {
      cosOrSin = &sin;
      cosOrSinPythonString = "math.sin";
    }
  }

  SlhaDsbHiggsVevFunctionoid::~SlhaDsbHiggsVevFunctionoid()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
