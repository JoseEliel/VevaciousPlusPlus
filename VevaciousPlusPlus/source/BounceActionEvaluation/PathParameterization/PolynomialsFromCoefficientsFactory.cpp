/*
 * PolynomialsFromCoefficientsFactory.cpp
 *
 *  Created on: Jul 21, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/PolynomialsFromCoefficientsFactory.hpp"

namespace VevaciousPlusPlus
{

  PolynomialsFromCoefficientsFactory::PolynomialsFromCoefficientsFactory(
                                                   size_t const numberOfFields,
                           size_t const numberOfVaryingCoefficientsPerField ) :
    TunnelPathFactory( numberOfFields,
                       std::vector< double >( ( numberOfFields
                                       * numberOfVaryingCoefficientsPerField ),
                                              0.0 ) ),
    falseVacuum( numberOfFields ),
    vacuumDifference( numberOfFields ),
    referenceFieldIndex( 0 ),
    coefficientsBeyondLinear( numberOfVaryingCoefficientsPerField )
  {
    // This constructor is just an initialization list.
  }

  PolynomialsFromCoefficientsFactory::~PolynomialsFromCoefficientsFactory()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
