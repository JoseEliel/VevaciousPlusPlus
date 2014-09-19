/*
 * TunnelPath.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/TunnelPath.hpp"

namespace VevaciousPlusPlus
{

  TunnelPath::TunnelPath( size_t const numberOfFields,
                          std::vector< double > const& pathParameterization,
                          double const pathTemperature ) :
    numberOfFields( numberOfFields ),
    pathParameterization( pathParameterization ),
    pathTemperature( pathTemperature ),
    nonZeroTemperature( pathTemperature > 0.0 )
  {
    // This constructor is just an initialization list.
  }

  TunnelPath::~TunnelPath()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
