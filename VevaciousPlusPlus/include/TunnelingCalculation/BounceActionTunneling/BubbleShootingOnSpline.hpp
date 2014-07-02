/*
 * BubbleShootingOnSpline.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BUBBLESHOOTINGONSPLINE_HPP_
#define BUBBLESHOOTINGONSPLINE_HPP_

#include "CommonIncludes.hpp"
#include "BounceActionCalculator.hpp"
#include "TunnelPath.hpp"

namespace VevaciousPlusPlus
{

  class BubbleShootingOnSpline : public BounceActionCalculator
  {
  public:
    BubbleShootingOnSpline( PotentialFunction const& potentialFunction,
                            std::string const& xmlArguments );
    virtual
    ~BubbleShootingOnSpline();


    // This sets radialStepSize and estimatedRadialMaximum based on the lengths
    // of the field vectors of the vacua and the tunneling scale from
    // potentialFunction. The smallest energy scale is inverted to give a
    // length which is used to set estimatedRadialMaximum, while radialStepSize
    // is given by lengthScaleResolution divided by the largest energy scale.
    virtual void ResetVacua( PotentialMinimum const& falseVacuum,
                             PotentialMinimum const& trueVacuum );

    // This sets up the bubble profile, numerically integrates the bounce
    // action over it, and then returns the bounce action along the path given
    // by tunnelPath. Either S_4, the dimensionless quantum bounce action
    // integrated over four dimensions, or S_3(T), the dimensionful (in GeV)
    // thermal bounce action integrated over three dimensions at temperature T,
    // should be returned: S_3(T) if the temperature T given by tunnelPath is
    // greater than 0.0, S_4 otherwise.
    virtual double operator()( TunnelPath const& tunnelPath ) const;

    // This plots the fields as functions of the bubble radial value in a file
    // called plotFilename in .eps format, with each field plotted in the color
    // given by fieldColors: the field with index i is plotted in the color
    // given by fieldColors[ i ]. An empty string indicates that the field
    // should not be plotted.
    virtual void PlotBounceConfiguration( TunnelPath const& tunnelPath,
                                 std::vector< std::string > const& fieldColors,
                                       std::string const& plotFilename ) const;


  protected:
    static double const radiusDifferenceThreshold;

    size_t numberOfPotentialSegments;
    double lengthScaleResolution;
    double radialStepSize;
    double estimatedRadialMaximum;
    size_t shootAttempts;
    double auxiliaryThreshold;


    // This returns a spline polynomial approximation of the potential along
    // the path given by tunnelPath. It also sets radialStepSize and
    // estimatedRadialMaximum based on the lengths of the field vectors of
    // the vacua on tunnelPath and the tunneling scale from potentialFunction.
    SplinePotential PotentialAlongPath( TunnelPath const& tunnelPath ) const;

    // This evaluates the bounce action density at the given point on the
    // bubble profile.
    double BounceActionDensity( SplinePotential const& potentialApproximation,
                                TunnelPath const& tunnelPath,
                      BubbleRadialValueDescription const& profilePoint ) const;
  };




  // This sets radialStepSize and estimatedRadialMaximum based on the lengths
  // of the field vectors of the vacua and the tunneling scale from
  // potentialFunction. The smallest energy scale is inverted to give a
  // length which is used to set estimatedRadialMaximum, while radialStepSize
  // is given by lengthScaleResolution divided by the largest energy scale.
  inline void
  BubbleShootingOnSpline::ResetVacua( PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    double currentEnergySquared(
                potentialFunction.ScaleSquaredRelevantToTunneling( falseVacuum,
                                                                trueVacuum ) );
    double largestEnergySquared( trueVacuum.LengthSquared() );
    double smallestEnergySquared( falseVacuum.LengthSquared() );
    if( largestEnergySquared < smallestEnergySquared )
    {
      std::swap( largestEnergySquared,
                 smallestEnergySquared );
    }
    if( smallestEnergySquared > currentEnergySquared )
    {
      smallestEnergySquared = currentEnergySquared;
    }
    else if( largestEnergySquared < currentEnergySquared )
    {
      largestEnergySquared = currentEnergySquared;
    }
    estimatedRadialMaximum = ( 1.0 / sqrt( smallestEnergySquared ) );
    radialStepSize = ( lengthScaleResolution / sqrt( largestEnergySquared ) );
  }

  // This evaluates the bounce action density at the given point on the
  // bubble profile.
  inline double BubbleShootingOnSpline::BounceActionDensity(
                                 SplinePotential const& potentialApproximation,
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
#endif /* BUBBLESHOOTINGONSPLINE_HPP_ */
