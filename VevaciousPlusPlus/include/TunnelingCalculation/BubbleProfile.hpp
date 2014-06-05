/*
 * BubbleProfile.hpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BUBBLEPROFILE_HPP_
#define BUBBLEPROFILE_HPP_

#include "../CommonIncludes.hpp"
#include "../PotentialEvaluation/SimplePolynomial.hpp"
#include "PathFieldsAndPotential.hpp"
#include "OdeintBubbleDerivatives.hpp"
#include "OdeintBubbleObserver.hpp"

namespace VevaciousPlusPlus
{

  class BubbleProfile
  {
  public:
    BubbleProfile( PathFieldsAndPotential const& pathFieldsAndPotential,
                   double const initialIntegrationStepSize,
                   double const initialIntegrationEndRadius );
    virtual
    ~BubbleProfile();


    // This tries the average of undershootAuxiliary and overshootAuxiliary
    // and uses it to update one of them based on whether its energy would
    // lead to undershooting or overshooting in the absence of the damping
    // term, and repeats this for a total of energyConservingUndershootAttempts
    // times.
    void UndampedUndershoot( size_t const energyConservingUndershootAttempts );

    // This tries to find the perfect shot undershootOvershootAttempts times,
    // then returns the bubble profile in terms of the auxiliary variable based
    // on the best shot. It integrates the auxiliary variable derivative to
    // increasing radial values until it is definite that the initial auxiliary
    // value gave an undershoot or an overshoot, or until the auxiliary value
    // at the largest radial value is within shootingThreshold of 0.
    std::vector< BubbleRadialValueDescription > const&
    DampedProfile( size_t const undershootOvershootAttempts,
                   double const shootingThreshold );


  protected:
    std::vector< BubbleRadialValueDescription > auxiliaryProfile;
    std::vector< BubbleRadialValueDescription > odeintProfile;
    PathFieldsAndPotential const& pathFieldsAndPotential;
    OdeintBubbleDerivatives bubbleDerivatives;
    OdeintBubbleObserver bubbleObserver;
    double integrationStepSize;
    double integrationEndRadius;
    double undershootAuxiliary;
    double overshootAuxiliary;
    double initialAuxiliary;
    std::vector< double > initialConditions;
    double const twiceDampingFactorPlusOne;
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
  };




  // This tries the average of undershootAuxiliary and overshootAuxiliary
  // and uses it to update one of them based on whether its energy would
  // lead to undershooting or overshooting in the absence of the damping
  // term, and repeats this for a total of energyConservingUndershootAttempts
  // times.
  inline void BubbleProfile::UndampedUndershoot(
                              size_t const energyConservingUndershootAttempts )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "BubbleProfile::UndampedUndershoot( "
    << energyConservingUndershootAttempts << " ) called:"
    << std::endl << "undershootAuxiliary = " << undershootAuxiliary
    << ", overshootAuxiliary = " << overshootAuxiliary;
    std::cout << std::endl;/**/
    for( size_t undershootGuessStep( 0 );
         undershootGuessStep < energyConservingUndershootAttempts;
         ++undershootGuessStep )
    {
      initialAuxiliary
      = ( 0.5 * ( undershootAuxiliary + overshootAuxiliary ) );
      if( pathFieldsAndPotential.PotentialApproximation( initialAuxiliary )
          < 0.0 )
      {
        overshootAuxiliary = initialAuxiliary;
      }
      else
      {
        undershootAuxiliary = initialAuxiliary;
      }
      // debugging:
      /**/std::cout << "tried currentAuxiliary = " << initialAuxiliary
      << ", V(p) = "
      << pathFieldsAndPotential.PotentialApproximation( initialAuxiliary )
      << ", now undershootAuxiliary = " << undershootAuxiliary
      << ", overshootAuxiliary = " << overshootAuxiliary;
      std::cout << std::endl;/**/
    }
    // At this point we reset overshootAuxiliary for the integration including
    // damping.
    overshootAuxiliary = 1.0;
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

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "BubbleProfile::ShootFromInitialConditions() called.";
    std::cout << std::endl;
    size_t integrationSteps =/**/
    boost::numeric::odeint::integrate( bubbleDerivatives,
                                       initialConditions,
                                       integrationStepSize,
                                       integrationEndRadius,
                                       integrationStepSize,
                                       bubbleObserver );

    RecordFromOdeintProfile();

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "integrationSteps = " << integrationSteps << ", bubble profile:";
    for( std::vector< BubbleRadialValueDescription >::const_iterator
         bubbleBit( auxiliaryProfile.begin() );
         bubbleBit < auxiliaryProfile.end();
         ++bubbleBit )
    {
      std::cout << "r = " << bubbleBit->radialValue << ", p = "
      << bubbleBit->auxiliaryValue << ", dp/dr = "
      << bubbleBit->auxiliarySlope << std::endl;
    }
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
#endif /* BUBBLEPROFILE_HPP_ */
