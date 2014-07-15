/*
 * TunnelPathFactory.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/TunnelPathFactory.hpp"

namespace VevaciousPlusPlus
{

  TunnelPathFactory::TunnelPathFactory( size_t const numberOfFields,
                          std::vector< double > const& zeroParameterization ) :
    numberOfFields( numberOfFields ),
    zeroParameterization( zeroParameterization )
  {
    // This constructor is just an initialization list.
  }

  TunnelPathFactory::~TunnelPathFactory()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
