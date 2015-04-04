/*
 * UndershootOvershootBubble.hpp
 *
 *  Created on: Nov 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef UNDERSHOOTOVERSHOOTBUBBLE_HPP_
#define UNDERSHOOTOVERSHOOTBUBBLE_HPP_

#include "CommonIncludes.hpp"
#include "boost/numeric/odeint/integrate/integrate.hpp"
#include "boost/math/special_functions/bessel.hpp"
#include "BasicFunctions/SimplePolynomial.hpp"
#include "SplinePotential.hpp"
#include "PathParameterization/TunnelPath.hpp"
#include "OdeintBubbleDerivatives.hpp"
#include "OdeintBubbleObserver.hpp"
#include "BubbleProfile.hpp"

namespace VevaciousPlusPlus
{

  class UndershootOvershootBubble : public BubbleProfile
  {
  public:
    UndershootOvershootBubble( double const initialIntegrationStepSize,
                               double const initialIntegrationEndRadius,
                               unsigned int const allowShootingAttempts,
                               double const shootingThreshold );
    virtual ~UndershootOvershootBubble();


    // This tries to find the perfect shot undershootOvershootAttempts times,
    // then sets auxiliaryProfile to be the bubble profile in terms of the
    // auxiliary variable based on the best shot. It integrates the auxiliary
    // variable derivative to increasing radial values until it is definite
    // that the initial auxiliary value gave an undershoot or an overshoot, or
    // until the auxiliary value at the largest radial value is within
    // shootingThreshold of auxiliaryAtRadialInfinity. It also correctly sets
    // auxiliaryAtRadialInfinity beforehand and afterwards
    // auxiliaryAtBubbleCenter.
    virtual void CalculateProfile( TunnelPath const& tunnelPath,
                                   SplinePotential const& pathPotential );

    // This interpolates between the 2 values of the auxiliary variable in
    // auxiliaryProfile with radial values either side of radialValue.
    virtual double AuxiliaryAt( double const radialValue ) const;

    // This interpolates between the 2 values of the auxiliary slope in
    // auxiliaryProfile with radial values either side of radialValue.
    virtual double AuxiliarySlopeAt( double const radialValue ) const;

    // This returns the maximum radial value up to which a profile can be
    // reasonably plotted.
    virtual double MaximumPlotRadius() const
    { return auxiliaryProfile.back().radialValue; }

    // This returns the discretized bubble profile calculated by the last call
    // of CalculateProfile.
    std::vector< BubbleRadialValueDescription > const& AuxiliaryProfile() const
    { return auxiliaryProfile; }

    // This just returns the path auxiliary at a radial value of 0.
    double AuxiliaryAtBubbleCenter() const{ return auxiliaryAtBubbleCenter; }

    // This just returns the path auxiliary at a radial value of infinity
    // (because of numerical effects, the path false vacuum might not be at
    // zero).
    double AuxiliaryAtRadialInfinity() const
    { return auxiliaryAtRadialInfinity; }


  protected:
    static double const auxiliaryPrecisionResolution;

    std::vector< BubbleRadialValueDescription > auxiliaryProfile;
    double auxiliaryAtBubbleCenter;
    double auxiliaryAtRadialInfinity;
    std::vector< BubbleRadialValueDescription > odeintProfile;
    double integrationStepSize;
    double integrationStartRadius;
    double integrationEndRadius;
    double undershootAuxiliary;
    double overshootAuxiliary;
    double initialAuxiliary;
    std::vector< double > initialConditions;
    double shootingThresholdSquared;
    unsigned int const allowShootingAttempts;
    bool worthIntegratingFurther;
    bool currentShotGoodEnough;
    TunnelPath const* tunnelPath;


    // This walks along auxiliaryProfile looking for the segment which starts
    // before radialValue and ends after it, in terms of the radial variable,
    // then returns the index of that segment and the difference between
    // radialValue and the radial start of the segment.
    std::pair< size_t, double >
    SegmentAndRemainder( double const radialValue ) const;

    // This returns the weights of the ends of the segment with index
    // segmentAndRemainder.first, according to how close
    // segmentAndRemainder.second is to each end, normalized to summing to 1.
    std::pair< double, double > SegmentEndWeights(
                std::pair< size_t, double > const& segmentAndRemainder ) const;

    // This looks through odeintProfile to see if there was a definite
    // undershoot or overshoot, setting undershootAuxiliary or
    // overshootAuxiliary respectively, as well as setting
    // worthIntegratingFurther. (It could sort odeintProfile based on radial
    // value, but odeint should have filled it in order.) Then it appends
    // odeintProfile to auxiliaryProfile.
    void RecordFromOdeintProfile( TunnelPath const& tunnelPath );

    // This performs the integration based on what is in initialConditions. It
    // also sets undershootAuxiliary, overshootAuxiliary, and
    // worthIntegratingFurther based on whether the integration showed that
    // there definitely was an undershoot, there definitely was an overshoot,
    // or that the auxiliary variable has not yet gotten within
    // shootingThresholdSquared^(1/2) of auxiliaryAtRadialInfinity.
    void ShootFromInitialConditions( TunnelPath const& tunnelPath,
                                     SplinePotential const& pathPotential );

    // This returns the slope of the solution for the bubble equation of motion
    // along the path in terms of p, which is either the derivative of
    // 2*sinh(x)/x for T != 0 or of 4*I_1(x)/x.
    double sinhOrBesselScaledSlope( bool const nonZeroTemperature,
                                    double const scaledRadius ) const;
  };




  // This interpolates between the 2 values of the auxiliary variable in
  // auxiliaryProfile with radial values either side of radialValue.
  inline double
  UndershootOvershootBubble::AuxiliaryAt( double const radialValue ) const
  {
    if( radialValue <= 0.0 )
    {
      return auxiliaryAtBubbleCenter;
    }
    else if( radialValue > auxiliaryProfile.back().radialValue )
    {
      return 0.0;
    }
    std::pair< size_t, double >
    segmentAndRemainder( SegmentAndRemainder( radialValue ) );
    std::pair< double, double >
    segmentEndWeights( SegmentEndWeights( segmentAndRemainder ) );
    return ( ( segmentEndWeights.first
               * auxiliaryProfile[ segmentAndRemainder.first ].auxiliaryValue )
             + ( segmentEndWeights.second
        * auxiliaryProfile[ segmentAndRemainder.first + 1 ].auxiliaryValue ) );
  }

  // This interpolates between the 2 values of the auxiliary slope in
  // auxiliaryProfile with radial values either side of radialValue.
  inline double
  UndershootOvershootBubble::AuxiliarySlopeAt( double const radialValue ) const
  {
    if( ( radialValue <= 0.0 )
        ||
        ( radialValue > auxiliaryProfile.back().radialValue ) )
    {
      return 0.0;
    }
    std::pair< size_t, double >
    segmentAndRemainder( SegmentAndRemainder( radialValue ) );
    std::pair< double, double >
    segmentEndWeights( SegmentEndWeights( segmentAndRemainder ) );
    return ( ( segmentEndWeights.first
               * auxiliaryProfile[ segmentAndRemainder.first ].auxiliarySlope )
             + ( segmentEndWeights.second
        * auxiliaryProfile[ segmentAndRemainder.first + 1 ].auxiliarySlope ) );
  }

  // This walks along auxiliaryProfile looking for the segment which starts
  // before radialValue and ends after it, in terms of the radial variable,
  // then returns the index of that segment and the difference between
  // radialValue and the radial start of the segment.
  inline std::pair< size_t, double >
  UndershootOvershootBubble::SegmentAndRemainder(
                                               double const radialValue ) const
  {
    size_t segmentIndex( 0 );
    size_t const maximumIndexForLoop( auxiliaryProfile.size() - 1 );
    while( ( segmentIndex < maximumIndexForLoop )
           &&
           ( radialValue > auxiliaryProfile[ ++segmentIndex ].radialValue ) ){}
    // The condition of the loop iterates segmentIndex until it is at the point
    // where auxiliaryProfile[ segmentIndex ] (after increment) starts beyond
    // radialValue, so we know that segmentIndex is 1 beyond what we want to
    // return.
    return std::pair< size_t, double >( segmentIndex - 1,
          ( radialValue - auxiliaryProfile[ segmentIndex - 1 ].radialValue ) );
  }

  // This returns the weights of the ends of the segment with index
  // segmentAndRemainder.first, according to how close
  // segmentAndRemainder.second is to each end, normalized to summing to 1.
  inline std::pair< double, double >
  UndershootOvershootBubble::SegmentEndWeights(
                 std::pair< size_t, double > const& segmentAndRemainder ) const
  {
    if( segmentAndRemainder.first >= ( auxiliaryProfile.size() - 1 ) )
    {
      return std::pair< double, double >( 1.0, 0.0 );
    }
    double const weightForLargerRadius( segmentAndRemainder.second /
                ( auxiliaryProfile[ segmentAndRemainder.first + 1 ].radialValue
               - auxiliaryProfile[ segmentAndRemainder.first ].radialValue ) );
    return std::pair< double, double >( ( 1.0 - weightForLargerRadius ),
                                        weightForLargerRadius );
  }

  // This performs the integration based on what is in initialConditions. It
  // also sets undershootAuxiliary, overshootAuxiliary, and
  // worthIntegratingFurther based on whether the integration showed that
  // there definitely was an undershoot, there definitely was an overshoot,
  // or that the auxiliary variable has not yet gotten within
  // shootingThresholdSquared^(1/2) of auxiliaryAtRadialInfinity.
  inline void UndershootOvershootBubble::ShootFromInitialConditions(
                                                  TunnelPath const& tunnelPath,
                                         SplinePotential const& pathPotential )
  {
    odeintProfile.clear();

    OdeintBubbleDerivatives bubbleDerivatives( pathPotential,
                                               tunnelPath );
    OdeintBubbleObserver bubbleObserver( odeintProfile );
    boost::numeric::odeint::integrate( bubbleDerivatives,
                                       initialConditions,
                                       integrationStartRadius,
                                       integrationEndRadius,
                                       integrationStepSize,
                                       bubbleObserver );
    RecordFromOdeintProfile( tunnelPath );
  }

  // This returns the slope of the solution for the bubble equation of motion
  // along the path in terms of p, which is either the derivative of
  // 2*sinh(x)/x for T != 0 or of 4*I_1(x)/x.
  inline double UndershootOvershootBubble::sinhOrBesselScaledSlope(
                                                 bool const nonZeroTemperature,
                                              double const scaledRadius ) const
  {
    if( nonZeroTemperature )
    {
      return ( ( 2.0 * ( ( scaledRadius * cosh( scaledRadius ) )
                         - sinh( scaledRadius ) ) )
               / ( scaledRadius * scaledRadius ) );
    }
    else
    {
      return ( ( 2.0 * ( ( scaledRadius *
                           ( boost::math::cyl_bessel_i( (int)0,
                                                        scaledRadius )
                             + boost::math::cyl_bessel_i( (int)2,
                                                          scaledRadius ) ) )
                         - ( 2.0 * boost::math::cyl_bessel_i( (int)1,
                                                           scaledRadius ) ) ) )
               / ( scaledRadius * scaledRadius ) );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* UNDERSHOOTOVERSHOOTBUBBLE_HPP_ */
