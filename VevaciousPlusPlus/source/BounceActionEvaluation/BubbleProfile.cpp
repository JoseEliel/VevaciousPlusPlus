/*
 * BubbleProfile.cpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BubbleProfile.hpp"

namespace VevaciousPlusPlus
{
  double const BubbleProfile::auxiliaryPrecisionResolution( 1.0e-6 );

  BubbleProfile::BubbleProfile( SplinePotential const& pathPotential,
                                TunnelPath const& tunnelPath,
                                double const initialIntegrationStepSize,
                                double const initialIntegrationEndRadius ) :
    auxiliaryProfile(),
    odeintProfile( 1,
                   BubbleRadialValueDescription() ),
    pathPotential( pathPotential ),
    tunnelPath( tunnelPath ),
    bubbleDerivatives( pathPotential,
                       tunnelPath ),
    bubbleObserver( odeintProfile ),
    integrationStepSize( initialIntegrationStepSize ),
    integrationStartRadius( initialIntegrationStepSize ),
    integrationEndRadius( initialIntegrationEndRadius ),
    undershootAuxiliary( NAN ),
    overshootAuxiliary( NAN ),
    initialAuxiliary( NAN ),
    initialConditions( 2 ),
    twoPlusTwiceDampingFactor( 8.0 ),
    shootingThresholdSquared( NAN ),
    shootAttemptsLeft( 0 ),
    worthIntegratingFurther( true ),
    currentShotGoodEnough( false )
  {
    if( tunnelPath.NonZeroTemperature() )
    {
      twoPlusTwiceDampingFactor = 6.0;
    }
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "BubbleProfile::BubbleProfile( pathPotential =" << std::endl;
    std::cout << pathPotential.AsDebuggingString() << std::endl;
    std::cout << "tunnelPath =" << std::endl;
    std::cout << tunnelPath.AsDebuggingString() << std::endl;
    std::cout << "initialIntegrationStepSize = " << initialIntegrationStepSize
    << std::endl;
    std::cout << "initialIntegrationEndRadius = "
    << initialIntegrationEndRadius << std::endl;
    std::cout << std::endl;*/
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
  BubbleProfile::operator()( size_t const undershootOvershootAttempts,
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
      undershootAuxiliary -= pathPotential.DefiniteOvershootAuxiliary();
    }
    if( overshootAuxiliary >= pathPotential.StartOfFinalSegment() )
    {
      overshootAuxiliary -= pathPotential.DefiniteOvershootAuxiliary();
    }

    // This loop is broken out of if the shoot attempt seems to have been close
    // enough that the integration would take too long to find an overshoot or
    // undershoot, or that the shot was dead on.
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "BubbleProfile::operator( " << undershootOvershootAttempts << ", "
    << shootingThreshold << " ) called. Initially, undershootAuxiliary = "
    << undershootAuxiliary << ", overshootAuxiliary = " << overshootAuxiliary;
    std::cout << std::endl;*/
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

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "undershootAuxiliary = " << undershootAuxiliary
      << ", overshootAuxiliary = " << overshootAuxiliary
      << " -> initialAuxiliary = " << initialAuxiliary;
      std::cout << std::endl;*/

      // We cannot start at r = 0, as the damping term is proportional to 1/r,
      // so the initial conditions are set by a Euler step assuming that near
      // r = 0, p goes as p_0 + p_2 r^2 (as the bubble should have smooth
      // fields at its center); hence d^2p/dr^2 (= 2 p_2) at r = 0 is
      // ( dV/dp ) / ( ( 1 + dampingFactor ) |df/dp|^2 ).
      // The initial step should be big enough that the initial conditions for
      // Boost::odeint will not suffer from precision problems from being too
      // close to the path panic minimum.
      if( initialAuxiliary <= 0.0 )
      {
        // If we're in the last segment, we go with a full solution of the
        // equation linearized in p:
        // [d^2/dr^2 + (dampingFactor/r)d/dr - ((d^2V/dp^2)/(|df/dp|^2))] p = 0
        // assuming that p is close enough to a minimum that dV/dp is
        // proportional to p and that dp/dr is small enough that we can neglect
        // the (dp/dr)^2 (df/dp).(d^2f/dp^2) part of (df/dp).(d^2f/dr^2).
        // The solutions are actually quite neat:
        // T = 0: p - p_0 = (p_i - p_0) * ( [ (4 * I_1(b*r)) / (b*r) ] - 1 )
        // T != 0: p - p_0 = (p_i - p_0) * ( [ (2 * sinh(b*r)) / (b*r) ] - 1 )
        // where I_1 is the modified Bessel function of the 1st kind,
        // b^2 = ((d^2V/dp^2)/(|df/dp|^2), p_0 is the auxiliary value at the
        // path panic minimum, and p_i is the initial value of the auxiliary
        // variable for the shoot.

        // We need the value of r such that |dp/dr| is larger than
        // auxiliaryPrecisionResolution but not much larger. We start assuming
        // that r is small and that we can expand out p(r) as p_0 + p_2 r^2,
        // which gives the same result as in the else statement complementary
        // to this branch, but here dV/dp = 2 * (p-p_0) * d^2V/dp^2.
        double const initialPositiveAuxiliary( initialAuxiliary
                                + pathPotential.DefiniteOvershootAuxiliary() );
        double const scaledSecondDerivative(
                pathPotential.SecondDerivativeNearPathPanic( initialAuxiliary )
                       / tunnelPath.SlopeSquared( initialPositiveAuxiliary ) );
        double const
        initialQuadraticCoefficient( 2.0 * initialAuxiliary
                                         * scaledSecondDerivative );

        // Because the potential is simply 2 minima with a maximum in between
        // (either originally so or truncated to it), and because
        // undershootAuxiliary is already past the maximum, the slope of the
        // potential at initialAuxiliary must be negative (in this case, a
        // negative initialAuxiliary times a positive second derivative), hence
        // the negative numerator.
        integrationStartRadius = ( -auxiliaryPrecisionResolution
               / ( 2.0 * initialQuadraticCoefficient * integrationStepSize ) );

        // debugging:
        /*std::cout << "in last segment. undershootAuxiliary = "
        << undershootAuxiliary << ", overshootAuxiliary = "
        << overshootAuxiliary << ", p = " << initialPositiveAuxiliary
        << " (" << pathPotential.DefiniteOvershootAuxiliary() << " - "
        << -initialAuxiliary << "), initialQuadraticCoefficient = "
        << initialQuadraticCoefficient << ", integrationStartRadius = "
        << integrationStartRadius;
        std::cout << std::endl;*/

        if( integrationStartRadius <= integrationStepSize )
        {
          // If integrationStartRadius turns out to be relatively small, we
          // carry on with the Euler step from the small-r approximation.
          initialConditions[ 0 ] = ( initialPositiveAuxiliary
                                     + ( initialQuadraticCoefficient
                                         * integrationStartRadius
                                         * integrationStartRadius ) );
          initialConditions[ 1 ] = ( 2.0 * initialQuadraticCoefficient
                                         * integrationStartRadius );

          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "small enough integrationStartRadius = " << integrationStartRadius
          << " <= integrationStepSize = " << integrationStepSize;
          std::cout << std::endl;*/
        }
        else
        {
          // If it turns out that maybe the initial step is too large to
          // consider the small-r expansion valid, we use the better
          // Bessel/sinh approximation mentioned above.
          double const inverseRadialScale( sqrt( scaledSecondDerivative ) );
          double const minimumScaledSlope( -auxiliaryPrecisionResolution
                                 / ( initialAuxiliary * inverseRadialScale ) );
          double scaledRadius( std::max( log( minimumScaledSlope ),
                              ( inverseRadialScale * integrationStepSize ) ) );
          double scaledSlope( sinhOrBesselScaledSlope(
                                               tunnelPath.NonZeroTemperature(),
                                                       scaledRadius ) );
          while( ( scaledSlope > ( 2.0 * minimumScaledSlope ) )
                 &&
                 ( scaledSlope > auxiliaryPrecisionResolution ) )
          {
            scaledRadius *= 0.5;
            scaledSlope
            = sinhOrBesselScaledSlope( tunnelPath.NonZeroTemperature(),
                                       scaledRadius );
            // debugging:
            /*std::cout << std::endl << "debugging:"
            << std::endl
            << "at scaled r = " << scaledRadius << ", scaled slope = "
            << scaledSlope;
            std::cout << std::endl;*/
          }
          while( scaledSlope < minimumScaledSlope )
          {
            scaledRadius = std::min( ( 2.0 * scaledRadius ),
                                     ( scaledRadius + 1.0 ) );
            scaledSlope
            = sinhOrBesselScaledSlope( tunnelPath.NonZeroTemperature(),
                                       scaledRadius );
            // debugging:
            /*std::cout << std::endl << "debugging:"
            << std::endl
            << "at scaled r = " << scaledRadius << ", scaled slope = "
            << scaledSlope;
            std::cout << std::endl;*/
          }
          // At this point, the slope of p(r) at
          // r = inverseRadialScale * scaledSlope ) should be large enough that
          // the magnitude of its product with integrationStepSize should be
          // larger than auxiliaryPrecisionResolution and thus the numeric
          // integration should be able to proceed normally.
          double sinhOrBesselPart( NAN );
          if( tunnelPath.NonZeroTemperature() )
          {
            sinhOrBesselPart = ( sinh( scaledRadius ) / scaledRadius );
          }
          else
          {
            sinhOrBesselPart = ( ( 2.0 * boost::math::cyl_bessel_i( (int)1,
                                                               scaledRadius ) )
                                 / scaledRadius );
          }
          initialConditions[ 0 ] = ( pathPotential.DefiniteOvershootAuxiliary()
                                     + ( initialAuxiliary
                                    * ( ( 2.0 * sinhOrBesselPart ) - 1.0 ) ) );
          initialConditions[ 1 ]
          = ( initialAuxiliary * inverseRadialScale * scaledSlope );
          integrationStartRadius = ( scaledRadius / inverseRadialScale );
          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "needed to do Bessel/sinh approximation: inverseRadialScale = "
          << inverseRadialScale << ", minimumScaledSlope = "
          << minimumScaledSlope << ", scaledRadius = " << scaledRadius
          << ", scaledSlope = " << scaledSlope << ", sinhOrBesselPart = "
          << sinhOrBesselPart;
          std::cout << std::endl;*/
        }
      }
      else
      {
        // If we're not starting in the last segment, we should not be
        // suffering from any of the problems due to having to roll very slowly
        // for a long r.
        double const initialPotentialDerivative( pathPotential.FirstDerivative(
                                                          initialAuxiliary ) );

        // debugging:
        /*std::cout << "not in last segment. undershootAuxiliary = "
        << undershootAuxiliary << ", overshootAuxiliary = "
        << overshootAuxiliary << ", trying p = " << initialAuxiliary << " ("
        << pathPotential.DefiniteOvershootAuxiliary() << " - "
        << ( pathPotential.DefiniteOvershootAuxiliary() - initialAuxiliary )
        << "), initialPotentialDerivative = " << initialPotentialDerivative;
        std::cout << std::endl;*/

        double const initialQuadraticCoefficient( initialPotentialDerivative
                                                  / ( twoPlusTwiceDampingFactor
                             * tunnelPath.SlopeSquared( initialAuxiliary ) ) );

        // Because the potential is simply 2 minima with a maximum in between
        // (either originally so or truncated to it), and because
        // undershootAuxiliary is already past the maximum, the slope of the
        // potential at initialAuxiliary must be negative, hence the negative
        // numerator.
        integrationStartRadius = ( -auxiliaryPrecisionResolution
               / ( 2.0 * initialQuadraticCoefficient * integrationStepSize ) );

        initialConditions[ 0 ] = ( initialAuxiliary
                                   + ( initialQuadraticCoefficient
                                       * integrationStartRadius
                                       * integrationStartRadius ) );
        initialConditions[ 1 ] = ( 2.0 * initialQuadraticCoefficient
                                       * integrationStartRadius );

        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "initialQuadraticCoefficient = " << initialQuadraticCoefficient
        << ", integrationStepSize = " << integrationStepSize
        << ", integrationStartRadius = " << integrationStartRadius;
        std::cout << std::endl;*/
      }


      // We have to ensure that the end radius is larger than the start radius.
      integrationEndRadius = std::max( integrationEndRadius,
                                       ( 2.0 * integrationStartRadius ) );


      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "after Euler 1st step, initialConditions = { "
      << initialConditions[ 0 ] << ", " << initialConditions[ 1 ]
      << " }. integrationStepSize = " << integrationStepSize
      << ", integrationStartRadius = " << integrationStartRadius
      << ", integrationEndRadius = " << integrationEndRadius
      << ", shootAttemptsLeft = " << shootAttemptsLeft;
      std::cout << std::endl;*/

      ShootFromInitialConditions();

      while( worthIntegratingFurther )
      {
        integrationStartRadius = auxiliaryProfile.back().radialValue;
        initialConditions[ 0 ] = auxiliaryProfile.back().auxiliaryValue;
        initialConditions[ 1 ] = auxiliaryProfile.back().auxiliarySlope;
        integrationEndRadius = ( 2.0 * integrationStartRadius );

        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "after start of loop integrating further, initialConditions = { "
        << initialConditions[ 0 ] << ", " << initialConditions[ 1 ]
        << " }. integrationStartRadius = " << integrationStartRadius
        << ", integrationEndRadius = " << integrationEndRadius;
        std::cout << std::endl;*/

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
      << ", " << tunnelPath.FieldsString( bubbleBit->auxiliaryValue )
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
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "overshoot! radialIndex = " << radialIndex
        << ", overshootAuxiliary = " << overshootAuxiliary << " -> "
        << initialAuxiliary;
        std::cout << std::endl;*/
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
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "undershoot! radialIndex = " << radialIndex
        << ", undershootAuxiliary = " << undershootAuxiliary << " -> "
        << initialAuxiliary;
        std::cout << std::endl;*/
        undershootAuxiliary = initialAuxiliary;
        worthIntegratingFurther = false;
        currentShotGoodEnough = false;
        break;
      }
      ++radialIndex;
    }
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "radialIndex = " << radialIndex;
    std::cout << std::endl;*/
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
      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "worthIntegratingFurther = " << worthIntegratingFurther
      << ", auxiliaryProfile.size() = " << auxiliaryProfile.size();
      std::cout << std::endl;*/
      std::vector< double > falseConfiguration( tunnelPath.NumberOfFields() );
      tunnelPath.PutOnPathAt( falseConfiguration,
                              0.0 );
      std::vector< double > currentConfiguration( falseConfiguration.size() );
      tunnelPath.PutOnPathAt( currentConfiguration,
                              auxiliaryProfile.back().auxiliaryValue );
      std::vector< double > initialConfiguration( falseConfiguration.size() );
      tunnelPath.PutOnPathAt( initialConfiguration,
                              initialAuxiliary );
      double initialDistanceSquared( 0.0 );
      double currentDistanceSquared( 0.0 );
      double fieldDifference( 0.0 );
      double falseVacuumField( 0.0 );
      for( size_t fieldIndex( 0 );
           fieldIndex < tunnelPath.NumberOfFields();
           ++fieldIndex )
      {
        falseVacuumField = falseConfiguration[ fieldIndex ];
        fieldDifference = ( currentConfiguration[ fieldIndex ]
                            - falseVacuumField );
        currentDistanceSquared += ( fieldDifference * fieldDifference );
        fieldDifference = ( initialConfiguration[ fieldIndex ]
                            - falseVacuumField );
        initialDistanceSquared += ( fieldDifference * fieldDifference );
      }
      currentShotGoodEnough
      = ( currentDistanceSquared
          < ( shootingThresholdSquared * initialDistanceSquared ) );
      worthIntegratingFurther = !currentShotGoodEnough;

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "currentDistanceSquared = " << currentDistanceSquared
      << ", initialDistanceSquared = " << initialDistanceSquared
      << ", shootingThresholdSquared = " << shootingThresholdSquared
      << ", currentShotGoodEnough = " << currentShotGoodEnough;
      std::cout << std::endl;*/
    }
  }

} /* namespace VevaciousPlusPlus */
