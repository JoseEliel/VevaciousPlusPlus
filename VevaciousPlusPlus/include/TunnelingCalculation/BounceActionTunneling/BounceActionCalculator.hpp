/*
 * BounceActionCalculator.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEACTIONCALCULATOR_HPP_
#define BOUNCEACTIONCALCULATOR_HPP_

#include "CommonIncludes.hpp"
#include "TunnelPath.hpp"

namespace VevaciousPlusPlus
{

  class BounceActionCalculator
  {
  public:
    BounceActionCalculator();
    virtual
    ~BounceActionCalculator();


    // This should calculate the bounce action along the path given by
    // tunnelPath.
    virtual double operator()( TunnelPath const& tunnelPath ) const = 0;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEACTIONCALCULATOR_HPP_ */
