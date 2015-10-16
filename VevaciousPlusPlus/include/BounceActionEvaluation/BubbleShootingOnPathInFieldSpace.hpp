/*
 * BubbleShootingOnPathInFieldSpace.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BUBBLESHOOTINGONPATHINFIELDSPACE_HPP_
#define BUBBLESHOOTINGONPATHINFIELDSPACE_HPP_

#include "CommonIncludes.hpp"
#include <limits>
#include "BounceActionCalculator.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PathParameterization/TunnelPath.hpp"
#include "OneDimensionalPotentialAlongPath.hpp"
#include "BubbleProfile.hpp"
#include "UndershootOvershootBubble.hpp"
#include "BubbleRadialValueDescription.hpp"

namespace VevaciousPlusPlus
{

  class BubbleShootingOnPathInFieldSpace : public BounceActionCalculator
  {
  public:
    BubbleShootingOnPathInFieldSpace(
                                    PotentialFunction const& potentialFunction,
                                      double const lengthScaleResolution,
                                      size_t const shootAttempts );
    virtual ~BubbleShootingOnPathInFieldSpace();


    // This sets radialStepSize to be lengthScaleResolution times the length
    // scale of the problem, which is taken to be 1.0 divided by the energy
    // scale. The energy scale is taken to be the square root of the scale
    // returned by potentialFunction.ScaleSquaredRelevantToTunneling( ... ) or
    // by tunnelingTemperature, whichever is larger. The maximum radius should
    // by infinity, but we cannot integrate to infinity, so we choose
    // estimatedRadialMaximum to be initially twice the length scale, as it
    // gets extended over the course of the calculation anyway if necessary.
    virtual void ResetVacua( PotentialMinimum const& falseVacuum,
                             PotentialMinimum const& trueVacuum,
                             double const tunnelingTemperature );

    // This sets up the bubble profile, numerically integrates the bounce
    // action over it, and then returns the bubble profile with calculated
    // bounce action along the path given by tunnelPath. Either S_4, the
    // dimensionless quantum bounce action integrated over four dimensions, or
    // S_3(T), the dimensionful (in GeV) thermal bounce action integrated over
    // three dimensions at temperature T, is calculated: S_3(T) if the
    // temperature T given by tunnelPath is greater than 0.0, S_4 otherwise.
    virtual BubbleProfile* operator()( TunnelPath const& tunnelPath,
                 OneDimensionalPotentialAlongPath const& pathPotential ) const;

    // This plots the fields as functions of the bubble radial value in a file
    // called plotFilename in .eps format, with each field plotted in the color
    // given by fieldColors: the field with index i is plotted in the color
    // given by fieldColors[ i ]. An empty string indicates that the field
    // should not be plotted.
    virtual void PlotBounceConfiguration( TunnelPath const& tunnelPath,
                                          BubbleProfile const& bubbleProfile,
                                 std::vector< std::string > const& fieldColors,
                                 unsigned int const plotResolution,
                                std::string const& gnuplotCommandIncludingPath,
                                       std::string const& plotFilename ) const;


  protected:
    static double const radiusDifferenceThreshold;

    double const lengthScaleResolution;
    double radialStepSize;
    double estimatedRadialMaximum;
    size_t const shootAttempts;
    double const auxiliaryThreshold;


    // This evaluates the bounce action density at the given point on the
    // bubble profile.
    double BounceActionDensity(
                OneDimensionalPotentialAlongPath const& potentialApproximation,
                                TunnelPath const& tunnelPath,
                      BubbleRadialValueDescription const& profilePoint ) const;
  };




  // This sets radialStepSize to be lengthScaleResolution times the length
  // scale of the problem, which is taken to be 1.0 divided by the energy
  // scale. The energy scale is taken to be the square root of the scale
  // returned by potentialFunction.ScaleSquaredRelevantToTunneling( ... ) or
  // by tunnelingTemperature, whichever is larger. The maximum radius should
  // by infinity, but we cannot integrate to infinity, so we choose
  // estimatedRadialMaximum to be initially twice the length scale, as it
  // gets extended over the course of the calculation anyway if necessary.
  inline void BubbleShootingOnPathInFieldSpace::ResetVacua(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                            double const tunnelingTemperature )
  {
    // We estimate the maximum radius to be twice the length scale, and the
    // step size to be lengthScaleResolution times the length scale.
    estimatedRadialMaximum = ( 2.0 / std::max( tunnelingTemperature,
          sqrt( potentialFunction.ScaleSquaredRelevantToTunneling( falseVacuum,
                                                            trueVacuum ) ) ) );
    radialStepSize = ( lengthScaleResolution * 0.5 * estimatedRadialMaximum );
  }

  // This evaluates the bounce action density at the given point on the
  // bubble profile.
  inline double BubbleShootingOnPathInFieldSpace::BounceActionDensity(
                OneDimensionalPotentialAlongPath const& potentialApproximation,
                                                  TunnelPath const& tunnelPath,
                       BubbleRadialValueDescription const& profilePoint ) const
  {
    double const currentAuxiliary( profilePoint.auxiliaryValue );
    double kineticTerm( profilePoint.auxiliarySlope );
    kineticTerm *= ( 0.5 * kineticTerm
                         * tunnelPath.SlopeSquared( currentAuxiliary ) );
    return ( kineticTerm + potentialApproximation( currentAuxiliary ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* BUBBLESHOOTINGONPATHINFIELDSPACE_HPP_ */