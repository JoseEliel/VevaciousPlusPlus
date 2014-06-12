/*
 * BubbleProfile.cpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{
  double const BubbleProfile::auxiliaryPrecisionResolution( 1.0e-6 );

  BubbleProfile::BubbleProfile(
                          PathFieldsAndPotential const& pathFieldsAndPotential,
                                double const initialIntegrationStepSize,
                                double const initialIntegrationEndRadius ) :
    auxiliaryProfile(),
    odeintProfile( 1,
                   BubbleRadialValueDescription() ),
    pathFieldsAndPotential( pathFieldsAndPotential ),
    pathPotential( pathFieldsAndPotential.PotentialApproximation() ),
    bubbleDerivatives( pathFieldsAndPotential ),
    bubbleObserver( odeintProfile ),
    integrationStepSize( initialIntegrationStepSize ),
    integrationStartRadius( initialIntegrationStepSize ),
    integrationEndRadius( initialIntegrationEndRadius ),
    undershootAuxiliary( NAN ),
    overshootAuxiliary( NAN ),
    initialAuxiliary( NAN ),
    initialConditions( 2 ),
    initialPositiveAuxiliary( NAN ),
    initialPotentialDerivative( NAN ),
    initialQuadraticCoefficient( NAN ),
    twiceDampingFactorPlusOne( ( 2.0 * pathFieldsAndPotential.DampingFactor() )
                               + 1.0 ),
    shootingThresholdSquared( NAN ),
    shootAttemptsLeft( 0 ),
    worthIntegratingFurther( true ),
    currentShotGoodEnough( false )
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
    shootAttemptsLeft = undershootOvershootAttempts;
    shootingThresholdSquared = ( shootingThreshold * shootingThreshold );
    undershootAuxiliary = pathPotential.DefiniteUndershootAuxiliary();
    overshootAuxiliary = pathPotential.DefiniteOvershootAuxiliary();

    // If undershootAuxiliary or overshootAuxiliary is in the final segment, it
    // is set to be the (negative) offset from the end of the final segment.
    // (This should mitigate precision issues for trying to start very close to
    // the path panic minimum.)
    if( undershootAuxiliary >= pathPotential.StartOfFinalSegment() )
    {
      undershootAuxiliary = ( pathPotential.DefiniteOvershootAuxiliary()
                              - undershootAuxiliary );
    }
    if( overshootAuxiliary >= pathPotential.StartOfFinalSegment() )
    {
      overshootAuxiliary = ( pathPotential.DefiniteOvershootAuxiliary()
                             - overshootAuxiliary );
    }

    // This loop is broken out of if the shoot attempt seems to have been close
    // enough that the integration would take too long to find an overshoot or
    // undershoot, or that the shot was dead on.
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "BubbleProfile::DampedProfile(" << undershootOvershootAttempts << ", "
    << shootingThreshold << " ) called:";
    std::cout << std::endl;/**/
    while( !currentShotGoodEnough
           &&
           ( shootAttemptsLeft > 0 ) )
    {
      worthIntegratingFurther = true;
      auxiliaryProfile.clear();
      integrationStartRadius = integrationStepSize;

      // It shouldn't ever happen that undershootAuxiliary is negative while
      // undershootAuxiliary is positive, as then the undershoot would be at a
      // larger auxiliary value than the overshoot.
      if( ( undershootAuxiliary > 0.0 )
          &&
          ( overshootAuxiliary <= 0.0 ) )
      {
        initialAuxiliary = ( 0.5 * ( undershootAuxiliary
                                     + overshootAuxiliary
                              + pathPotential.DefiniteOvershootAuxiliary() ) );
        if( initialAuxiliary >= pathPotential.StartOfFinalSegment() )
        {
          initialAuxiliary = ( pathPotential.StartOfFinalSegment()
                              - initialAuxiliary );
        }
      }
      else
      {
        initialAuxiliary
        = ( 0.5 * ( undershootAuxiliary + overshootAuxiliary ) );
      }

      // We cannot start at r = 0, as the damping term is proportional to 1/r,
      // so the initial conditions are set by a Euler step assuming that near
      // r = 0, p goes as p_0 + p_2 r^2 (as the bubble should have smooth
      // fields at its center); hence d^2p/dr^2 at r = 0 is
      // ( dV/dp ) / ( ( 1 + 2 * dampingFactor ) |df/dp|^2 ).
      // The initial step should be big enough that the initial conditions for
      // Boost::odeint will not suffer from precision problems from being too
      // close to the path panic minimum.
      if( initialAuxiliary <= 0.0 )
      {
        initialPositiveAuxiliary = ( pathPotential.DefiniteOvershootAuxiliary()
                                     + initialAuxiliary );
        initialPotentialDerivative
        = pathPotential.DerivativeNearPathPanic( initialAuxiliary );
        // debugging:
        /**/std::cout << "undershootAuxiliary = " << undershootAuxiliary
        << ", overshootAuxiliary = " << overshootAuxiliary << ", trying p = "
        << initialPositiveAuxiliary << " ("
        << pathPotential.DefiniteOvershootAuxiliary() << " - "
        << -initialAuxiliary << "), initialPotentialDerivative = "
        << initialPotentialDerivative;
        std::cout << std::endl;/**/
      }
      else
      {
        initialPositiveAuxiliary = initialAuxiliary;
        initialPotentialDerivative
        = pathPotential.FirstDerivative( initialAuxiliary );
        // debugging:
        /**/std::cout << "undershootAuxiliary = " << undershootAuxiliary
        << ", overshootAuxiliary = " << overshootAuxiliary << ", trying p = "
        << initialAuxiliary << " ("
        << pathPotential.DefiniteOvershootAuxiliary() << " - "
        << ( pathPotential.DefiniteOvershootAuxiliary() - initialAuxiliary )
        << "), initialPotentialDerivative = " << initialPotentialDerivative;
        std::cout << std::endl;/**/
      }
      initialQuadraticCoefficient = ( initialPotentialDerivative
                                      / ( twiceDampingFactorPlusOne
                              * pathFieldsAndPotential.FieldDerivativesSquared(
                                                initialPositiveAuxiliary ) ) );

      integrationStartRadius = ( -auxiliaryPrecisionResolution
                     / ( initialQuadraticCoefficient * integrationStepSize ) );


      initialConditions[ 0 ] = ( initialPositiveAuxiliary
                                 + ( initialQuadraticCoefficient
                                     * integrationStartRadius
                                     * integrationStartRadius ) );
      initialConditions[ 1 ] = ( 2.0 * initialQuadraticCoefficient
                                     * integrationStartRadius );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "after Euler 1st step, initialConditions = { "
      << initialConditions[ 0 ] << ", " << initialConditions[ 1 ]
      << " }. integrationStepSize = " << integrationStepSize
      << ", integrationStartRadius = " << integrationStartRadius
      << ", integrationEndRadius = " << integrationEndRadius;
      std::cout << std::endl;/**/

      ShootFromInitialConditions();

      while( worthIntegratingFurther )
      {
        integrationStartRadius = auxiliaryProfile.back().radialValue;
        initialConditions[ 0 ] = auxiliaryProfile.back().auxiliaryValue;
        initialConditions[ 1 ] = auxiliaryProfile.back().auxiliarySlope;
        integrationEndRadius = ( integrationStartRadius
                                 + integrationStartRadius );

        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "after start of loop integrating further, initialConditions = { "
        << initialConditions[ 0 ] << ", " << initialConditions[ 1 ]
        << " }. integrationStartRadius = " << integrationStartRadius
        << ", integrationEndRadius = " << integrationEndRadius;
        std::cout << std::endl;/**/

        ShootFromInitialConditions();
      }
      --shootAttemptsLeft;
    }
    // At the end of the loop, initialAuxiliary is either within
    // 2^(-undershootOvershootAttempts) of p_crit, or was close enough that the
    // integration to decide if it was an undershoot or overshoot would take
    // too long.

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Final bubble profile:" << std::endl;
    for( std::vector< BubbleRadialValueDescription >::const_iterator
         bubbleBit( auxiliaryProfile.begin() );
         bubbleBit < auxiliaryProfile.end();
         ++bubbleBit )
    {
      std::cout << "r = " << bubbleBit->radialValue << ", p = "
      << bubbleBit->auxiliaryValue << ", dp/dr = " << bubbleBit->auxiliarySlope
      << ", "
      << pathFieldsAndPotential.FieldsString( bubbleBit->auxiliarySlope )
      << std::endl;
    }
    std::cout << std::endl;/**/

    return auxiliaryProfile;
  }

  // This looks through odeintProfile to see if there was a definite
  // undershoot or overshoot, setting undershootAuxiliary or
  // overshootAuxiliary respectively, as well as setting
  // worthIntegratingFurther. (It could sort odeintProfile based on radial
  // value, but odeint should have filled it in order.) Then it appends
  // odeintProfile to auxiliaryProfile.
  void BubbleProfile::RecordFromOdeintProfile()
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "BubbleProfile::RecordFromOdeintProfile() called."
    << " odeintProfile.size() = " << odeintProfile.size();
    std::cout << std::endl;/**/
    // We start from the beginning of odeintProfile so that we record only as
    // much of the bubble profile as there is before the shot starts to roll
    // backwards or overshoot.
    size_t radialIndex( 0 );
    while( radialIndex < odeintProfile.size() )
    {
      // If the shot has gone past the false vacuum, it was definitely an
      // overshoot.
      if( odeintProfile[ radialIndex ].auxiliaryValue < 0.0 )
      {
        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "overshoot! radialIndex = " << radialIndex
        << ", overshootAuxiliary = " << overshootAuxiliary << " -> "
        << initialAuxiliary;
        std::cout << std::endl;/**/
        overshootAuxiliary = initialAuxiliary;
        worthIntegratingFurther = false;
        currentShotGoodEnough = false;
        break;
      }
      // If the shot is rolling backwards without having yet reached the false
      // vacuum, it was definitely an undershoot.
      else if( odeintProfile[ radialIndex ].auxiliarySlope > 0.0 )
      {
        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "undershoot! radialIndex = " << radialIndex
        << ", undershootAuxiliary = " << undershootAuxiliary << " -> "
        << initialAuxiliary;
        std::cout << std::endl;/**/
        undershootAuxiliary = initialAuxiliary;
        worthIntegratingFurther = false;
        currentShotGoodEnough = false;
        break;
      }
      ++radialIndex;
    }
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "radialIndex = " << radialIndex;
    std::cout << std::endl;/**/
    if( radialIndex < odeintProfile.size() )
    {
      auxiliaryProfile.insert( auxiliaryProfile.end(),
                               ( odeintProfile.begin() + 1 ),
                               ( odeintProfile.begin() + radialIndex ) );
    }
    else
    {
      auxiliaryProfile.insert( auxiliaryProfile.end(),
                               ( odeintProfile.begin() + 1 ),
                               odeintProfile.end() );
    }
    odeintProfile.clear();
    // If there wasn't an undershoot or overshoot, currentShotGoodEnough
    // has to be set based on whether the shot got close enough to the false
    // vacuum.
    if( worthIntegratingFurther )
    {
      double initialDistanceSquared( 0.0 );
      double currentDistanceSquared( 0.0 );
      double fieldDifference( 0.0 );
      double falseVacuumField( 0.0 );
      double const currentAuxiliary( auxiliaryProfile.back().auxiliaryValue );
      for( std::vector< SimplePolynomial >::const_iterator
           fieldValue( pathFieldsAndPotential.FieldPath().begin() );
           fieldValue < pathFieldsAndPotential.FieldPath().end();
           ++fieldValue )
      {
        falseVacuumField = fieldValue->CoefficientVector().front();
        fieldDifference = ( (*fieldValue)( currentAuxiliary )
                            - falseVacuumField );
        currentDistanceSquared += ( fieldDifference * fieldDifference );
        fieldDifference = ( (*fieldValue)( initialAuxiliary )
                            - falseVacuumField );
        initialDistanceSquared += ( fieldDifference * fieldDifference );
      }
      currentShotGoodEnough
      = ( currentDistanceSquared
          < ( shootingThresholdSquared * initialDistanceSquared ) );
      worthIntegratingFurther = !currentShotGoodEnough;

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "currentDistanceSquared = " << currentDistanceSquared
      << ", initialDistanceSquared = " << initialDistanceSquared
      << ", shootingThresholdSquared = " << shootingThresholdSquared
      << ", currentShotGoodEnough = " << currentShotGoodEnough;
      std::cout << std::endl;/**/
    }
  }

} /* namespace VevaciousPlusPlus */
