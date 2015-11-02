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
                                             size_t const indexInValuesVector,
                     SlhaSourcedParameterFunctionoid const& vevEuclideanLength,
                                SlhaSourcedParameterFunctionoid const& tanBeta,
                                                       bool const sinNotCos ) :
    SlhaSourcedParameterFunctionoid( indexInValuesVector ),
    vevEuclideanLength( vevEuclideanLength ),
    vevIndex( vevEuclideanLength.IndexInValuesVector() ),
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
