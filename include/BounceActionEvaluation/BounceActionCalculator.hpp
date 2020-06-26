/*
 * BounceActionCalculator.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEACTIONCALCULATOR_HPP_
#define BOUNCEACTIONCALCULATOR_HPP_

#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "BubbleProfile.hpp"
#include "PathParameterization/TunnelPath.hpp"
#include "OneDimensionalPotentialAlongPath.hpp"

namespace VevaciousPlusPlus
{

  class BounceActionCalculator
  {
  public:
    BounceActionCalculator() {}

    virtual ~BounceActionCalculator() {}


    // This should prepare the BounceActionCalculator for a path or set of
    // paths between the given vacua, possibly because it might need to set up
    // things based on characteristic energy scales for example. By default it
    // does nothing.
    virtual void ResetVacua( PotentialFunction const& potentialFunction,
                             PotentialMinimum const& falseVacuum,
                             PotentialMinimum const& trueVacuum,
                             double const tunnelingTemperature ) {}

    // This should calculate the bubble profile and bounce action along the
    // path given by tunnelPath and return them in a new BubbleProfile object,
    // which must have its memory managed by the calling code (ideally a
    // unique_ptr< BubbleProfile > would be returned, but we're sticking to
    // allowing non-C++11-compliant compilers). Either S_4, the dimensionless
    // quantum bounce action integrated over four dimensions, or S_3(T), the
    // dimensionful (in GeV) thermal bounce action integrated over three
    // dimensions at temperature T, should be set in the returned
    // BubbleProfile: S_3(T) if the temperature T given by tunnelPath is
    // greater than 0.0, S_4 otherwise.
    virtual BubbleProfile* operator()( TunnelPath const& tunnelPath,
             OneDimensionalPotentialAlongPath const& pathPotential ) const = 0;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEACTIONCALCULATOR_HPP_ */
