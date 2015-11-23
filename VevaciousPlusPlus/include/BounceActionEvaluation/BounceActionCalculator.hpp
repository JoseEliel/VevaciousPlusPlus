/*
 * BounceActionCalculator.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEACTIONCALCULATOR_HPP_
#define BOUNCEACTIONCALCULATOR_HPP_

#include "CommonIncludes.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PathParameterization/TunnelPath.hpp"
#include "OneDimensionalPotentialAlongPath.hpp"
#include "BubbleProfile.hpp"

namespace VevaciousPlusPlus
{

  class BounceActionCalculator
  {
  public:
    BounceActionCalculator();
    virtual ~BounceActionCalculator();


    // This should prepare the BounceActionCalculator for a path or set of
    // paths between the given vacua, possibly because it might need to set up
    // things based on characteristic energy scales for example. By default it
    // does nothing.
    virtual void ResetVacua( PotentialFunction const& potentialFunction,
                             PotentialMinimum const& falseVacuum,
                             PotentialMinimum const& trueVacuum,
                             double const tunnelingTemperature ){}

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

    // This plots the fields as functions of the spatial variables in a file
    // called plotFilename in .eps format, with each field plotted in the color
    // given by fieldColors: the field with index i is plotted in the color
    // given by fieldColors[ i ]. An empty string indicates that the field
    // should not be plotted.
    virtual void PlotBounceConfiguration( TunnelPath const& tunnelPath,
                                          BubbleProfile const& bubbleProfile,
                                  std::vector< std::string > const& fieldNames,
                                 std::vector< std::string > const& fieldColors,
                                          std::string const& plotFilename,
                                     unsigned int const plotResolution ) const;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEACTIONCALCULATOR_HPP_ */
