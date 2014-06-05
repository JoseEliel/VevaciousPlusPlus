/*
 * BubbleProfile.cpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{

  BubbleProfile::BubbleProfile(
                          PathFieldsAndPotential const& pathFieldsAndPotential,
                                double const initialIntegrationStepSize,
                                double const initialIntegrationEndRadius ) :
    auxiliaryProfile(),
    odeintProfile(),
    pathFieldsAndPotential( pathFieldsAndPotential ),
    bubbleDerivatives( pathFieldsAndPotential ),
    bubbleObserver( odeintProfile ),
    integrationStepSize( initialIntegrationStepSize ),
    integrationEndRadius( initialIntegrationEndRadius ),
    undershootAuxiliary( 0.0 ),
    overshootAuxiliary( 1.0 ),
    currentAuxiliary( NAN ),
    initialConditions( 2 ),
    twiceDampingFactorPlusOne( ( 2.0 * pathFieldsAndPotential.DampingFactor() )
                               + 1.0 ),
    shootingThresholdSquared( NAN ),
    worthIntegratingFurther( false )
  {
    // This constructor is just an initialization list.
  }

  BubbleProfile::~BubbleProfile()
  {
    // This does nothing.
  }


  // This tries to find the perfect shot undershootOvershootAttempts times,
  // then returns the bubble profile in terms of the auxiliary variable based
  // on the best shot. It integrates the auxiliary variable derivative to
  // increasing radial values until it is definite that the initial auxiliary
  // value gave an undershoot or an overshoot, or until the auxiliary value
  // at the largest radial value is within shootingThreshold of 0.
  std::vector< BubbleRadialValueDescription > const&
  BubbleProfile::DampedProfile( size_t const undershootOvershootAttempts,
                                double const shootingThreshold )
  {
    shootingThresholdSquared = ( shootingThreshold * shootingThreshold );
    size_t shootAttempts( 0 );

    // This loop is broken out of if the shoot attempt seems to have been close
    // enough that the integration would take too long to find an overshoot or
    // undershoot, or that the shot was dead on.
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "BubbleProfile::DampedProfile(" << undershootOvershootAttempts << ", "
    << shootingThreshold << " ) called:";
    std::cout << std::endl;/**/
    while( StillSearching( ++shootAttempts ) )
    {
      auxiliaryProfile.clear();
      currentAuxiliary
      = ( 0.5 * ( undershootAuxiliary + overshootAuxiliary ) );
      // debugging:
      /**/std::cout << "undershootAuxiliary = " << undershootAuxiliary
      << ", overshootAuxiliary = " << overshootAuxiliary << ", trying p = "
      << currentAuxiliary;
      std::cout << std::endl;/**/

      // We cannot start at r = 0, as the damping term is proportional to 1/r,
      // so the initial conditions are set by a Euler step assuming that near
      // r = 0, p goes as p_0 + p_2 r^2 (as the bubble should have smooth
      // fields at its center); hence d^2p/dr^2 at r = 0 is
      // ( dV/dp ) / ( ( 1 + 2 * dampingFactor ) |df/dp|^2 ).
      initialConditions[ 0 ] = currentAuxiliary;
      initialConditions[ 1 ]
      = ( ( bubbleDerivatives.PotentialDerivative( currentAuxiliary )
            * integrationStepSize )
          / ( twiceDampingFactorPlusOne
              * pathFieldsAndPotential.FieldDerivativesSquared(
                                                        currentAuxiliary ) ) );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "after Euler 1st step, initialConditions = { "
      << initialConditions[ 0 ] << ", " << initialConditions[ 1 ]
      << " }. integrationStepSize = " << integrationStepSize
      << ", integrationEndRadius = " << integrationEndRadius;
      std::cout << std::endl;/**/

      ShootFromInitialConditions();


      while( worthIntegratingFurther )
      {
        initialConditions[ 0 ] = auxiliaryProfile.back().radialValue;
        initialConditions[ 1 ] = auxiliaryProfile.back().auxiliaryValue;

        integrationEndRadius = ( 2.0 * initialConditions[ 0 ] );
        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "after start of loop integrating further, initialConditions = { "
        << initialConditions[ 0 ] << ", " << initialConditions[ 1 ]
        << " }. integrationEndRadius = " << integrationEndRadius;
        std::cout << std::endl;/**/

        ShootFromInitialConditions();
      }
      // Now we have to decide if initialConditions[ 0 ] is the new
      // undershootAuxiliary or overshootAuxiliary, or if we need to extend
      // integrationRadius, or if the current guess overshoots with a
      // negligible amount of kinetic energy.

      if( odeintBubbleObserver.DefinitelyUndershot() )
      {
        undershootAuxiliary = currentAuxiliary;
      }
      else if( odeintBubbleObserver.DefinitelyOvershot() )
      {
        overshootAuxiliary = currentAuxiliary;
      }
      else
      {
        // We break out of the loop leaving currentAuxiliary as it is if the
        // shot was good enough.
        break;
      }
      currentAuxiliary
      = ( 0.5 * ( undershootAuxiliary + overshootAuxiliary ) );
    }
    // At the end of the loop, currentAuxiliary is either within
    // 2^(-undershootOvershootAttempts) of p_crit, or was close enough that the
    // integration to decide if it was an undershot or overshot would take too
    // long.

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Final bubble profile:"
    << std::endl << odeintBubbleObserver.AsDebuggingString();
    std::cout << std::endl;/**/


    return auxiliaryProfile;
  }

} /* namespace VevaciousPlusPlus */
