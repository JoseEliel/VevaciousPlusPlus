/*
 * GradientMinimizer.hpp
 *
 *  Created on: Jun 30, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef GRADIENTMINIMIZER_HPP_
#define GRADIENTMINIMIZER_HPP_

#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimum.hpp"
#include <vector>

namespace VevaciousPlusPlus
{

  class GradientMinimizer
  {
  public:
    GradientMinimizer( PotentialFunction const& potentialFunction ) :
      potentialFunction( potentialFunction ) {}

    virtual ~GradientMinimizer() {}


    // This should find the minimum using startingPoint as a starting point,
    // obviously.
    virtual PotentialMinimum
    operator()( std::vector< double > const& startingPoint ) const = 0;

    // This should ensure that the minimizations are calculated at the given
    // temperature.
    virtual void SetTemperature( double const minimizationTemperature ) = 0;


  protected:
    PotentialFunction const& potentialFunction;
  };

} /* namespace VevaciousPlusPlus */
#endif /* GRADIENTMINIMIZER_HPP_ */
