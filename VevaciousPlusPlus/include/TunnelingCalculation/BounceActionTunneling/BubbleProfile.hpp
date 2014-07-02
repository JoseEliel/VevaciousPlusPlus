/*
 * BubbleProfile.hpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BUBBLEPROFILE_HPP_
#define BUBBLEPROFILE_HPP_

#include "CommonIncludes.hpp"
#include "boost/numeric/odeint/integrate/integrate.hpp"
#include "boost/math/special_functions/bessel.hpp"
#include "BasicFunctions/SimplePolynomial.hpp"
#include "SplinePotential.hpp"
#include "TunnelPath.hpp"
#include "OdeintBubbleDerivatives.hpp"
#include "OdeintBubbleObserver.hpp"

namespace VevaciousPlusPlus
{

  class BubbleProfile
  {
  public:
    BubbleProfile( SplinePotential const& potentialApproximation,
                   TunnelPath const& tunnelPath,
                   double const initialIntegrationStepSize,
                   double const initialIntegrationEndRadius );
    virtual
    ~BubbleProfile();


    // This tries to find the perfect shot undershootOvershootAttempts times,
    // then returns the bubble profile in terms of the auxiliary variable based
    // on the best shot. It integrates the auxiliary variable derivative to
    // increasing radial values until it is definite that the initial auxiliary
    // value gave an undershoot or an overshoot, or until the auxiliary value
    // at the largest radial value is within shootingThreshold of 0.
    std::vector< BubbleRadialValueDescription > const&
    operator()( size_t const undershootOvershootAttempts,
                double const shootingThreshold );

    // This returns the value that the auxiliary variable should have at the
    // center of the bubble.
    double AuxiliaryAtBubbleCenter() const;


  protected:
    static double const auxiliaryPrecisionResolution;

    std::vector< BubbleRadialValueDescription > auxiliaryProfile;
    std::vector< BubbleRadialValueDescription > odeintProfile;
    SplinePotential const& pathPotential;
    TunnelPath const& tunnelPath;
    OdeintBubbleDerivatives bubbleDerivatives;
    OdeintBubbleObserver bubbleObserver;
    double integrationStepSize;
    double integrationStartRadius;
    double integrationEndRadius;
    double undershootAuxiliary;
    double overshootAuxiliary;
    double initialAuxiliary;
    std::vector< double > initialConditions;
    double twoPlusTwiceDampingFactor;
    double shootingThresholdSquared;
    size_t shootAttemptsLeft;
    bool worthIntegratingFurther;
    bool currentShotGoodEnough;


    // This looks through odeintProfile to see if there was a definite
    // undershoot or overshoot, setting undershootAuxiliary or
    // overshootAuxiliary respectively, as well as setting
    // worthIntegratingFurther. (It could sort odeintProfile based on radial
    // value, but odeint should have filled it in order.) Then it appends
    // odeintProfile to auxiliaryProfile.
    void RecordFromOdeintProfile();

    // This performs the integration based on what is in initialConditions. It
    // also sets undershootAuxiliary, overshootAuxiliary, and
    // worthIntegratingFurther based on whether the integration showed that
    // there definitely was an undershoot, there definitely was an overshoot,
    // or that the auxiliary variable has not yet gotten within
    // shootingThresholdSquared^(1/2) of 0.
    void ShootFromInitialConditions();

    // This returns the slope of the solution for the bubble equation of motion
    // along the path in terms of p, which is either the derivative of
    // 2*sinh(x)/x for T != 0 or of 4*I_1(x)/x.
    double sinhOrBesselScaledSlope( bool const nonZeroTemperature,
                                    double const scaledRadius ) const;
  };





  // This returns the value that the auxiliary variable should have at the
  // center of the bubble.
  inline double BubbleProfile::AuxiliaryAtBubbleCenter() const
  {
    if( initialAuxiliary < 0.0 )
    {
      return ( pathPotential.DefiniteOvershootAuxiliary() + initialAuxiliary );
    }
    else
    {
      return initialAuxiliary;
    }
  }

  // This performs the integration based on what is in initialConditions. It
  // also sets undershootAuxiliary, overshootAuxiliary, and
  // worthIntegratingFurther based on whether the integration showed that
  // there definitely was an undershoot, there definitely was an overshoot,
  // or that the auxiliary variable has not yet gotten within
  // shootingThresholdSquared^(1/2) of 0.
  inline void BubbleProfile::ShootFromInitialConditions()
  {
    odeintProfile.clear();
    boost::numeric::odeint::integrate( bubbleDerivatives,
                                       initialConditions,
                                       integrationStartRadius,
                                       integrationEndRadius,
                                       integrationStepSize,
                                       bubbleObserver );

    RecordFromOdeintProfile();
  }

  // This returns the slope of the solution for the bubble equation of motion
  // along the path in terms of p, which is either the derivative of
  // 2*sinh(x)/x for T != 0 or of 4*I_1(x)/x.
  inline double
  BubbleProfile::sinhOrBesselScaledSlope( bool const nonZeroTemperature,
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
#endif /* BUBBLEPROFILE_HPP_ */
