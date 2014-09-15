/*
 * BubbleShootingOnSpline.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BubbleShootingOnSpline.hpp"

namespace VevaciousPlusPlus
{
  double const BubbleShootingOnSpline::radiusDifferenceThreshold( 0.01 );

  BubbleShootingOnSpline::BubbleShootingOnSpline(
                                    PotentialFunction const& potentialFunction,
                                        size_t const numberOfPotentialSegments,
                                            double const lengthScaleResolution,
                                                 size_t const shootAttempts ) :
    BounceActionCalculator( potentialFunction ),
    numberOfPotentialSegments( numberOfPotentialSegments ),
    lengthScaleResolution( lengthScaleResolution ),
    radialStepSize( NAN ),
    estimatedRadialMaximum( NAN ),
    shootAttempts( shootAttempts ),
    auxiliaryThreshold( 1.0E-6 )
  {
    // This constructor is just an initialization list.
  }

  BubbleShootingOnSpline::~BubbleShootingOnSpline()
  {
    // This does nothing.
  }


  // This sets up the bubble profile, numerically integrates the bounce
  // action over it, and then returns the bounce action along the path given
  // by tunnelPath. Either S_4, the dimensionless quantum bounce action
  // integrated over four dimensions, or S_3(T), the dimensionful (in GeV)
  // thermal bounce action integrated over three dimensions at temperature T,
  // should be returned: S_3(T) if the temperature T given by tunnelPath is
  // greater than 0.0, S_4 otherwise.
  double
  BubbleShootingOnSpline::operator()( TunnelPath const& tunnelPath ) const
  {
    SplinePotential potentialApproximation( PotentialAlongPath( tunnelPath ) );
    bool const nonZeroTemperature( tunnelPath.NonZeroTemperature() );
    BubbleProfile bubbleProfile( potentialApproximation,
                                 tunnelPath,
                                 radialStepSize,
                                 estimatedRadialMaximum );
    std::vector< BubbleRadialValueDescription > const&
    auxiliaryProfile( bubbleProfile( shootAttempts,
                                     auxiliaryThreshold ) );

    // We have a set of radial values r_i, path auxiliary values p(r_i), and
    // slopes dp/dr|_{r=r_i}, and can easily evaluate a set of "bounce action
    // densities" B_i = B(r_i). The numerical integral is then the sum of
    // B_i * [differential volume at r_i], which normally would be
    // B_i * r_i^dampingFactor * [solid angle] * ( r_{i+1} - r_i ).
    // In an attempt to reduce the numerical error, one can instead take
    // 0.5 * (B_i + B_{i+1}) * r_i^dampingFactor
    //     * [solid angle] * ( r_{i+1} - r_i ), or other averages that are
    // formally identical as the difference in radius vanishes. We take the
    // average B(r) over each differential volume:
    // ( 0.5 * (B_i + B_{i+1}) ) * [differential volume]
    // where [differential volume] is
    // ( r_i^dampingFactor * [solid angle] * ( r_{i+1} - r_i ) )
    // if ( r_{i+1} - r_i ) is small compared to r_{i+1} and r_i, or we take
    // [solid angle] * ( r_{i+1}^d - r_i^d )/d where d = (dampingFactor+1) if
    // the difference in radius is not small enough.
    // Since then say B_4 is counted with 2 volumes (with a factor of 0.5 each
    // time), we actually sum up each B_i with the sum of its 2 adjacent
    // volumes, hence we sum up
    // 0.5 * [solid angle] * B_i * r_i^(d-1) * ( r_{i+1} - r_{i-1} )
    // or 0.5 * [solid angle] * B_i * ( r_{i+1}^d - r_{i-1}^d )/d.
    // (Actually, we leave the common factor of 0.5 * [solid angle] until after
    // the end of the loop.)
    // The first and last contributions have to be treated slightly
    // differently, because of the lack of radii at indices beyond the range of
    // auxiliaryProfile.

    double previousRadius( 0.0 );
    double previousVolume( 0.0 );
    double currentRadius( auxiliaryProfile.front().radialValue );
    double currentVolume( currentRadius * currentRadius * currentRadius );
    double nextRadius( auxiliaryProfile[ 1 ].radialValue );
    double nextVolume( nextRadius * nextRadius * nextRadius );
    if( nonZeroTemperature )
    {
      currentVolume /= 3.0;
      nextVolume /= 3.0;
    }
    else
    {
      currentVolume *= ( 0.25 * currentRadius );
      nextVolume *= ( 0.25 * nextRadius );
    }
    // The bounce action up to the radius of auxiliaryProfile[ 1 ] is given by
    // B_{-1} (the bounce action density at r = 0.0, which has no kinetic
    // contribution because the bubble is smooth at its center by construction)
    // and B_0 as
    // ( 0.5 * ( B_{-1} + B_{0} ) * [volume from r = 0.0 to r_0] )
    // + ( 0.5 * ( B_{0} + B_{1} ) * [volume from r = r_0 to r_1] ).
    // The contribution from B_{1} is taken care of with the 1st iteration of
    // the loop, so the contribution up to the loop is given by
    // ( 0.5 * [solid angle] * B_{-1} * currentVolume )
    // + ( 0.5 * [solid angle] * B_{0} * nextVolume )
    // with the values currently in currentVolume and nextVolume, also with the
    // common factor of 0.5 * [solid angle] being left until after the loop.
    double
    bounceAction( ( currentVolume
                    * potentialApproximation(
                                    bubbleProfile.AuxiliaryAtBubbleCenter() ) )
                  + ( nextVolume
                      * ( BounceActionDensity( potentialApproximation,
                                               tunnelPath,
                                              auxiliaryProfile.front() ) ) ) );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "Before loop: currentRadius = " << currentRadius << ", currentVolume = "
    << currentVolume << ", nextRadius = " << nextRadius << ", nextVolume = "
    << nextVolume << ", p(0) = " << bubbleProfile.AuxiliaryAtBubbleCenter()
    << ", B_{-1} = "
    << potentialApproximation( bubbleProfile.AuxiliaryAtBubbleCenter() )
    << ", B_{0} = " << BounceActionDensity( potentialApproximation,
                                            tunnelPath,
                                            auxiliaryProfile.front() )
    << ", bounceAction = " << bounceAction;
    std::cout << std::endl;*/

    for( size_t radiusIndex( 1 );
         radiusIndex < ( auxiliaryProfile.size() - 1 );
         ++radiusIndex )
    {
      previousRadius = currentRadius;
      previousVolume = currentVolume;
      currentRadius = nextRadius;
      currentVolume = nextVolume;
      nextRadius = auxiliaryProfile[ radiusIndex + 1 ].radialValue;
      nextVolume = ( nextRadius * nextRadius * nextRadius );
      if( nonZeroTemperature )
      {
        nextVolume /= 3.0;
      }
      else
      {
        nextVolume *= ( 0.25 * nextRadius );
      }

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "previousRadius = " << previousRadius
      << ", previousVolume = " << previousVolume
      << ", currentRadius = " << currentRadius
      << ", currentVolume = " << currentVolume
      << ", nextRadius = " << nextRadius
      << ", nextVolume = " << nextVolume
      << ", B_" << radiusIndex << " = "
      << BounceActionDensity( potentialApproximation,
                              tunnelPath,
                              auxiliaryProfile[ radiusIndex ] );
      std::cout << std::endl;*/

      // B_i * r_i^(d-1) * ( r_{i+1} - r_{i-1} )
      // or B_i * ( r_{i+1}^d - r_{i-1}^d )/d.
      if( ( nextRadius - previousRadius )
          > ( radiusDifferenceThreshold * currentRadius ) )
      {
        bounceAction
        += ( BounceActionDensity( potentialApproximation,
                                  tunnelPath,
                                  auxiliaryProfile[ radiusIndex ] )
            * ( nextVolume - previousVolume ) );
      }
      else
      {
        double currentArea( currentRadius * currentRadius );
        if( !nonZeroTemperature )
        {
          currentArea *= currentRadius;
        }
        bounceAction
        += ( BounceActionDensity( potentialApproximation,
                                  tunnelPath,
                                  auxiliaryProfile[ radiusIndex ] )
             * ( nextRadius - previousRadius )
             * currentArea );
      }
      // The common factor of 0.5 * [solid angle] is being left until after the
      // loop.
      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "bounceAction = " << bounceAction;
      std::cout << std::endl;*/
    }
    // Now we add the last shell:
    double const currentAuxiliary( auxiliaryProfile.back().auxiliaryValue );
    double kineticTerm( auxiliaryProfile.back().auxiliarySlope );
    kineticTerm *= ( 0.5 * kineticTerm
                         * tunnelPath.SlopeSquared( currentAuxiliary ) );
    double const potentialTerm( potentialApproximation( currentAuxiliary ) );
    if( ( nextRadius - currentRadius )
        > ( radiusDifferenceThreshold * currentRadius ) )
    {
      bounceAction
      += ( ( kineticTerm + potentialTerm ) * ( nextVolume - currentVolume ) );
    }
    else
    {
      double currentArea( currentRadius * currentRadius );
      if( !nonZeroTemperature )
      {
        currentArea *= currentRadius;
      }
      bounceAction += ( ( kineticTerm + potentialTerm )
                        * ( nextRadius - currentRadius ) * currentArea );
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "After last shell, bounceAction = " << bounceAction;
    std::cout << std::endl;*/

    // Near the false vacuum at p = 0, the potential should be of the form
    // constant + p^2 * (d^2V/dp^2) / 2, so the bubble equations of motion can
    // be linearized to
    // [d^2/dr^2 + ([damping factor]/r) d/dr - ((d^2V/dp^2)/(|df/dp|^2))] p = 0
    // assuming that dp/dr is small enough that we can neglect the
    // (dp/dr)^2 (df/dp).(d^2f/dp^2) part of (df/dp).(d^2f/dr^2).
    // The solution for p for the quantum case of [damping factor] = 3 which
    // tends to 0 as r tends to infinity is proportional to ( K[1,(b*r)] / r ),
    // where K is the modified Bessel function of the 2nd kind and
    // b^2 = (d^2V/dp^2)/(|df/dp|^2)).
    // The solution for p for the thermal case of [damping factor] = 2 which
    // tends to 0 as r tends to infinity is ( exp[b(R-r)] R p(R) / r ).
    // The integrals of the potential term from R to infinity both in the
    // quantum case and in the thermal case can be done, yielding that the
    // potential integral is equal to the value of the potential term at R
    // times a factor (R^3 / b) * ( K[2, (b R)] / K[1, (b R)] ) for the quantum
    // case or (R / b^2 ) * ( 1 + b R ) for the thermal case.
    // The thermal kinetic term can also be integrated, yielding the value of
    // the kinetic term at R times a factor
    // ( R^3 * [ ( 1 + [(b R)/2] ) / ( 1 + b R )^2 ] ).
    // Unfortunately the quantum kinetic term cannot be integrated in a closed
    // form (in the sense that we have easy access to a numeric value), so we
    // consider the case where R is so large that the damping term can be
    // neglected, and get that d^2f/dr^2 = f d^2V/df^2 for a field f, which is
    // solved by f being proportional to exp(-r). Therefore the kinetic term
    // decays like e^( 2 - ((2 r)/R) ). The integral from R to infinity is
    // r^3 e(-(2r)/R) -> (19/8) R^4 / e^2, which combines with the e^2 from
    // e^( 2 - ((2 r)/R) ) to give the factor as simply 2.375 R^4.
    // (All solutions found with Mathematica 8.)

    // This is b in the mathematics above.
    double const
    inverseScale( sqrt( potentialApproximation.SecondDerivative( 0.0 ) ) );
    double const scaledRadius( inverseScale * nextRadius );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "inverseScale = " << inverseScale << ", scaledRadius = "
    << scaledRadius;
    std::cout << std::endl;*/

    if( nonZeroTemperature )
    {
      // Potential term: (R / b^2 ) * ( 1 + b R ).
      // Kinetic term: ( R^3 * [ ( 1 + [(b R)/2] ) / ( 1 + b R )^2 ] ).
      double const onePlusRatio( 1.0 + scaledRadius );
      bounceAction
      += ( nextRadius
           * ( ( ( potentialTerm * onePlusRatio )
                 / ( inverseScale * inverseScale ) )
               + ( ( kineticTerm * nextRadius * nextRadius
                     * ( 1.0 + ( 0.5 * scaledRadius ) ) )
                   / ( onePlusRatio * onePlusRatio ) ) ) );
    }
    else
    {
      // Potential term: (R^3 / b) * ( K[2, (b R)] / K[1, (b R)] ).
      // Kinetic term: (19/8) R^4.
      // The ratio of the Bessel functions tends to 1 quite quickly (for values
      // of bR > 100, the ratio is within 1.5% of 1) while the actual Bessel
      // functions drop very quickly and can easily have overflow errors in the
      // exponent.
      double potentialTimesFactor( potentialTerm / inverseScale );
      if( scaledRadius < 100.0 )
      {
        potentialTimesFactor *= ( boost::math::cyl_bessel_k( (int)2,
                                                             scaledRadius )
                                 / boost::math::cyl_bessel_k( (int)1,
                                                              scaledRadius ) );
      }
      bounceAction += ( nextRadius * nextRadius * nextRadius
                        * ( potentialTimesFactor
                            + ( 2.375 * kineticTerm * nextRadius ) ) );
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "After Bessel or exponent functions to infinity, bounceAction = "
    << bounceAction;
    std::cout << std::endl;*/

    // The common factor of 1/2 is combined with the solid angle of
    // 2 pi^2 (quantum) or 4 pi (thermal):
    if( nonZeroTemperature )
    {
      return ( bounceAction * 2.0 * boost::math::double_constants::pi );
    }
    else
    {
      return ( bounceAction * boost::math::double_constants::pi
                            * boost::math::double_constants::pi );
    }
  }

  // This plots the fields as functions of the bubble radial value in a file
  // called plotFilename in .eps format, with each field plotted in the color
  // given by fieldColors: the field with index i is plotted in the color
  // given by fieldColors[ i ]. An empty string indicates that the field
  // should not be plotted.
  void BubbleShootingOnSpline::PlotBounceConfiguration(
                                                  TunnelPath const& tunnelPath,
                                 std::vector< std::string > const& fieldColors,
                                        std::string const& plotFilename ) const
  {
    SplinePotential potentialApproximation( PotentialAlongPath( tunnelPath ) );
    BubbleProfile bubbleProfile( potentialApproximation,
                                 tunnelPath,
                                 radialStepSize,
                                 estimatedRadialMaximum );
    std::vector< BubbleRadialValueDescription > const&
    auxiliaryProfile( bubbleProfile( shootAttempts,
                                     auxiliaryThreshold ) );
    BOL::TwoDimensionalDataPlotter bubblePlotter( "/opt/local/bin/gnuplot",
                                                  plotFilename );

    std::cout << std::endl
    << "Path parameterized as:" << std::endl
    << tunnelPath.AsDebuggingString();
    std::cout << std::endl
    << "Potential parameterized as:" << std::endl
    << potentialApproximation.AsDebuggingString();
    std::cout << std::endl;
    std::cout << "Bubble profile:" << std::endl;
    for( std::vector< BubbleRadialValueDescription >::const_iterator
         bubbleBit( auxiliaryProfile.begin() );
         bubbleBit < auxiliaryProfile.end();
         ++bubbleBit )
    {
      std::cout << "r = " << bubbleBit->radialValue << ", p = "
      << bubbleBit->auxiliaryValue << ", dp/dr = " << bubbleBit->auxiliarySlope
      << ", "
      << tunnelPath.FieldsString( bubbleBit->auxiliaryValue )
      << std::endl;
    }
    std::cout << std::endl;

    size_t numberOfPlottedFields( std::min( fieldColors.size(),
                                potentialFunction.NumberOfFieldVariables() ) );
    BOL::TwoDimensionalDataPlotter::PlotDataVector plotData;
    BOL::TwoDimensionalDataPlotter::DoublePairVectorWithStringPair fieldData;
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfPlottedFields;
         ++fieldIndex )
    {
      if( !(fieldColors[ fieldIndex ].empty()) )
      {
        fieldData.second.first.assign( fieldColors[ fieldIndex ] );
        fieldData.second.second.assign(
                                potentialFunction.FieldNames()[ fieldIndex ] );
        plotData.push_back( fieldData );
      }
    }
    numberOfPlottedFields = plotData.size();
    std::vector< double >
    fieldConfiguration( potentialFunction.NumberOfFieldVariables() );
    for( size_t radiusIndex( 0 );
         radiusIndex < auxiliaryProfile.size();
         ++radiusIndex )
    {
      double const radialValue( auxiliaryProfile[ radiusIndex ].radialValue );
      tunnelPath.PutOnPathAt( fieldConfiguration,
                              auxiliaryProfile[ radiusIndex ].auxiliaryValue );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfPlottedFields;
           ++fieldIndex )
      {
        plotData[ fieldIndex ].first.push_back( std::make_pair( radialValue,
                                          fieldConfiguration[ fieldIndex ] ) );
      }
    }
    bubblePlotter.plotData( plotData,
                            "rho/(1/GeV)",
                            "field/GeV" );
  }

  // This returns a spline polynomial approximation of the potential along
  // the path given by tunnelPath.
  SplinePotential BubbleShootingOnSpline::PotentialAlongPath(
                                           TunnelPath const& tunnelPath ) const
  {
    SplinePotential potentialApproximation;
    std::vector< double >
    fieldConfiguration( potentialFunction.NumberOfFieldVariables() );
    tunnelPath.PutOnPathAt( fieldConfiguration,
                            0.0 );
    double const falseVacuumPotential( potentialFunction( fieldConfiguration,
                                             tunnelPath.TemperatureValue() ) );
    // We choose to take the cos of a linear distribution of equally-spaced
    // points so that we have a finer resolution of the potential near the
    // minima.
    double auxiliaryValue( 0.0 );
    for( size_t splinePoint( 1 );
         splinePoint < numberOfPotentialSegments;
         ++splinePoint )
    {
      auxiliaryValue = ( 0.5 * ( 1.0 - cos( ( splinePoint
                                          * boost::math::double_constants::pi )
                    / static_cast< double >( numberOfPotentialSegments ) ) ) );
      // The above might run into numerical precision issues if too many spline
      // segments are asked for, but it's probably OK up to even 10^5 segments
      // (cos(0.0) - cos(10^(-5)) = 5E-11 which should still be well-resolved
      // by doubles).
      // auxiliaryValue = ( static_cast< double >( splinePoint )
      //               / static_cast< double >( numberOfSplinesInPotential ) );
      tunnelPath.PutOnPathAt( fieldConfiguration,
                              auxiliaryValue );

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "auxiliaryValue = " << auxiliaryValue << ", fieldConfiguration = "
      << potentialFunction.FieldConfigurationAsMathematica(
                                                          fieldConfiguration );
      std::cout << std::endl;*/

      potentialApproximation.AddPoint( auxiliaryValue,
                                       ( potentialFunction( fieldConfiguration,
                                                tunnelPath.TemperatureValue() )
                                         - falseVacuumPotential ) );
    }
    tunnelPath.PutOnPathAt( fieldConfiguration,
                            1.0 );
    potentialApproximation.SetSpline( potentialFunction( fieldConfiguration,
                                                tunnelPath.TemperatureValue() )
                                      - falseVacuumPotential );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "BubbleShootingOnSpline::PotentialAlongPath( tunnelPath =" << std::endl
    << tunnelPath.AsDebuggingString() << " ) returning" << std::endl
    << potentialApproximation.AsDebuggingString();
    std::cout << std::endl;/**/
    return potentialApproximation;
  }

} /* namespace VevaciousPlusPlus */
