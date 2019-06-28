/*
 * UndershootOvershootBubble.cpp
 *
 *  Created on: Nov 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/UndershootOvershootBubble.hpp"

namespace VevaciousPlusPlus
{
  double const
  UndershootOvershootBubble::auxiliaryPrecisionResolution( 1.0e-7 );

  UndershootOvershootBubble::UndershootOvershootBubble(
                                       double const initialIntegrationStepSize,
                                      double const initialIntegrationEndRadius,
                                      unsigned int const allowShootingAttempts,
                                             double const shootingThreshold ) :
    BubbleProfile(),
    auxiliaryProfile(),
    auxiliaryAtBubbleCenter( -1.0 ),
    auxiliaryAtRadialInfinity( -1.0 ),
    odeintProfile( 1,
                   BubbleRadialValueDescription() ),
    integrationStepSize( initialIntegrationStepSize ),
    integrationStartRadius( initialIntegrationStepSize ),
    integrationEndRadius( initialIntegrationEndRadius ),
    undershootAuxiliary( 0.0 ),
    overshootAuxiliary( 1.0 ),
    initialAuxiliary( 0.5 ),
    initialConditions( 2 ),
    shootingThresholdSquared( shootingThreshold * shootingThreshold ),
    allowShootingAttempts( allowShootingAttempts ),
    worthIntegratingFurther( true ),
    currentShotGoodEnough( false ),
    badInitialConditions( false ),
    tunnelPath( NULL )
  {
    // This constructor is just an initialization list.
  }

  UndershootOvershootBubble::~UndershootOvershootBubble()
  {
    // This does nothing.
  }


  // This tries to find the perfect shot undershootOvershootAttempts times,
  // then returns the bubble profile in terms of the auxiliary variable based
  // on the best shot. It integrates the auxiliary variable derivative to
  // increasing radial values until it is definite that the initial auxiliary
  // value gave an undershoot or an overshoot, or until the auxiliary value
  // at the largest radial value is within shootingThreshold of
  // auxiliaryAtRadialInfinity. It also correctly sets
  // auxiliaryAtRadialInfinity beforehand and afterwards
  // auxiliaryAtBubbleCenter.
  void
  UndershootOvershootBubble::CalculateProfile( TunnelPath const& tunnelPath,
                        OneDimensionalPotentialAlongPath const& pathPotential )
  {
    auxiliaryAtRadialInfinity = pathPotential.AuxiliaryOfPathFalseVacuum();
    double twoPlusTwiceDampingFactor( 8.0 );
    if( tunnelPath.NonZeroTemperature() )
    {
      twoPlusTwiceDampingFactor = 6.0;
    }

    unsigned int shootAttemptsLeft( allowShootingAttempts );

    undershootAuxiliary = pathPotential.DefiniteUndershootAuxiliary();
    overshootAuxiliary = pathPotential.AuxiliaryOfPathPanicVacuum();

//    std::cout << "Just starting " << "Undershoot: " << undershootAuxiliary << " " << "Overshoot: " << overshootAuxiliary<< std::endl;

    // If undershootAuxiliary is too close to the path panic
    // vacuum, it is set to be the (negative) offset from the panic vacuum.
    // (This should mitigate precision issues for trying to start very close to
    // the path panic minimum.)
//    std::cout << "Threshold " << pathPotential.ThresholdForNearPathPanic() << std::endl;

    if( undershootAuxiliary >= pathPotential.ThresholdForNearPathPanic() )
    {
        undershootAuxiliary -= pathPotential.AuxiliaryOfPathPanicVacuum();
    }


    // This loop is broken out of if the shoot attempt seems to have been close
    // enough that the integration would take too long to find an overshoot or
    // undershoot, or that the shot was dead on.
    while( !currentShotGoodEnough
           &&
           ( shootAttemptsLeft > 0 ) )
    {
      worthIntegratingFurther = true;
      auxiliaryProfile.clear();
      integrationStartRadius = integrationStepSize;

//      std::cout << "integrationStartRadius I " << integrationStartRadius << std::endl;
//
//
//      std::cout << " " << "Undershoot: " << undershootAuxiliary << " " << "Overshoot: " << overshootAuxiliary<< std::endl;

      // It shouldn't ever happen that undershootAuxiliary is negative while
      // overshootAuxiliary is positive, as then the undershoot would be at a
      // larger auxiliary value than the overshoot.
      if( ( undershootAuxiliary > 0.0 )
          &&
          ( overshootAuxiliary <= 0.0 ) )
      {
        initialAuxiliary = ( 0.5 * ( undershootAuxiliary
                                     + overshootAuxiliary
                              + pathPotential.AuxiliaryOfPathPanicVacuum() ) );
        if( initialAuxiliary >= pathPotential.ThresholdForNearPathPanic() )
        {
          initialAuxiliary -= pathPotential.AuxiliaryOfPathPanicVacuum();
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
      // fields at its center); hence d^2p/dr^2 (= 2 p_2) at r = 0 is
      // ( dV/dp ) / ( ( 1 + dampingFactor ) |df/dp|^2 ).
      // The initial step should be big enough that the initial conditions for
      // Boost::odeint will not suffer from precision problems from being too
      // close to the path panic minimum.
      if( initialAuxiliary <= 0.0 )
      {
        // If we're close to the path panic minimum, we go with a full solution
        // of the equation linearized in p:
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
                                + pathPotential.AuxiliaryOfPathPanicVacuum() );
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

//        std::cout << "integrationStartRadius II " << integrationStartRadius << std::endl;

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
          }
          while( scaledSlope < minimumScaledSlope )
          {
            scaledRadius = std::min( ( 2.0 * scaledRadius ),
                                     ( scaledRadius + 1.0 ) );
            scaledSlope
            = sinhOrBesselScaledSlope( tunnelPath.NonZeroTemperature(),
                                       scaledRadius );
          }
          // At this point, the slope of p(r) at
          // r = ( inverseRadialScale * scaledSlope ) should be large enough
          // that the magnitude of its product with integrationStepSize should
          // be larger than auxiliaryPrecisionResolution and thus the numeric
          // integration should be able to proceed normally.
          double sinhOrBesselPart( -2.0 );
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
          // As above:
          // p - p_0 = (p_i - p_0) * ( [ (4 * I_1(b*r)) / (b*r) ] - 1 )
          // or p - p_0 = (p_i - p_0) * ( [ (2 * sinh(b*r)) / (b*r) ] - 1 )
          initialConditions[ 0 ] = ( pathPotential.AuxiliaryOfPathPanicVacuum()
                                     + ( initialAuxiliary
                                    * ( ( 2.0 * sinhOrBesselPart ) - 1.0 ) ) );
          initialConditions[ 1 ]
          = ( initialAuxiliary * inverseRadialScale * scaledSlope );
          integrationStartRadius = ( scaledRadius / inverseRadialScale );
//          std::cout << "integrationStartRadius III " << integrationStartRadius << std::endl;

        }
      }
      else
      {
        // If we're not starting in the last segment, we should not be
        // suffering from any of the problems due to having to roll very slowly
        // for a long r.
        double const initialPotentialDerivative( pathPotential.FirstDerivative(
                                                          initialAuxiliary ) );
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

//        std::cout << "integrationStartRadius IV " << integrationStartRadius << std::endl;

          initialConditions[ 0 ] = ( initialAuxiliary
                                   + ( initialQuadraticCoefficient
                                       * integrationStartRadius
                                       * integrationStartRadius ) );
        initialConditions[ 1 ] = ( 2.0 * initialQuadraticCoefficient
                                       * integrationStartRadius );
      }

      // We have to ensure that the end radius is larger than the start radius.
      integrationEndRadius = std::max( integrationEndRadius,
                                       ( 2.0 * integrationStartRadius ) );
      ShootFromInitialConditions( tunnelPath,
                                  pathPotential );

      while( worthIntegratingFurther )
      {
        integrationStartRadius = auxiliaryProfile.back().radialValue;
        initialConditions[ 0 ] = auxiliaryProfile.back().auxiliaryValue;
        initialConditions[ 1 ] = auxiliaryProfile.back().auxiliarySlope;
        integrationEndRadius = ( 2.0 * integrationStartRadius );

//        std::cout << "Initial Conditions V " << initialConditions[0] << " " << initialConditions[1] << " " << auxiliaryProfile.back().auxiliaryValue << " " << auxiliaryAtRadialInfinity << std::endl;
//
//        std::cout << "integrationStartRadius V " << integrationStartRadius << std::endl;

        ShootFromInitialConditions( tunnelPath,
                                    pathPotential );
      }
      --shootAttemptsLeft;
    }
    // At the end of the loop, initialAuxiliary is either within
    // 2^(-undershootOvershootAttempts) of p_crit, or was close enough that the
    // integration to decide if it was an undershoot or overshoot would take
    // too long.

    if( initialAuxiliary < 0.0 )
    {
      auxiliaryAtBubbleCenter = ( pathPotential.AuxiliaryOfPathPanicVacuum()
                                  + initialAuxiliary );
    }
    else
    {
      auxiliaryAtBubbleCenter = initialAuxiliary;
    }
  }

  // This looks through odeintProfile to see if there was a definite
  // undershoot or overshoot, setting undershootAuxiliary or
  // overshootAuxiliary respectively, as well as setting
  // worthIntegratingFurther. (It could sort odeintProfile based on radial
  // value, but odeint should have filled it in order.) Then it appends
  // odeintProfile to auxiliaryProfile.
  void UndershootOvershootBubble::RecordFromOdeintProfile(
                                                 TunnelPath const& tunnelPath )
  {
    // We start from the beginning of odeintProfile so that we record only as
    // much of the bubble profile as there is before the shot starts to roll
    // backwards or overshoot.
    size_t radialIndex( 0 );
    while( radialIndex < odeintProfile.size() )
    {
      // If the shot has gone past the false vacuum, it was definitely an
      // overshoot.
      if( odeintProfile[ radialIndex ].auxiliaryValue
          < auxiliaryAtRadialInfinity )
      {
//          std::cout<< "we are at overshoot detection "<< odeintProfile[ radialIndex ].auxiliaryValue <<" "<< auxiliaryAtRadialInfinity  <<std::endl;
//          std::cout << "throwing point " << initialAuxiliary <<std::endl;
//          std::cout << "Previous Undershoot " << undershootAuxiliary <<std::endl;
//          std::cout << "Previous Overshoot " << overshootAuxiliary <<std::endl;
        overshootAuxiliary = initialAuxiliary;
        worthIntegratingFurther = false;
        currentShotGoodEnough = false;
        break;
      }
      // If the shot is rolling backwards without having yet reached the false
      // vacuum, it was definitely an undershoot.
      else if( odeintProfile[ radialIndex ].auxiliarySlope > 0.0 )
      {
//        std::cout<< "we are at undershoot detection "<< odeintProfile[ radialIndex ].auxiliarySlope<< std::endl;
//        std::cout << "throwing point " << initialAuxiliary <<std::endl;
//        std::cout << "Previous Undershoot " << undershootAuxiliary <<std::endl;
//        std::cout << "Previous Overshoot " << overshootAuxiliary <<std::endl;
        undershootAuxiliary = initialAuxiliary;
        worthIntegratingFurther = false;
        currentShotGoodEnough = false;
        break;
      }
      ++radialIndex;
    }

    if( radialIndex < odeintProfile.size() )
    {
      if( radialIndex > 1 )
      {
        badInitialConditions = false;
        auxiliaryProfile.insert( auxiliaryProfile.end(),
                                 ( odeintProfile.begin() + 1 ),
                                 ( odeintProfile.begin() + radialIndex ) );
      }
      else
      {
          // Here we check whether the initial radius was so large that we end up at step 0 in a definite
          // undershoot/overshoot. In that case, we go back and set the initial step radius to be smaller.
          // this happens in ShootFromInitialConditions.
          badInitialConditions = true;
          std::cout<< " Rescaling initial integration radius in under/overshoot to help with"
                       << " detected numerical problems. Shooting again now."<<std::endl;

      }
    }
    else
    {
      badInitialConditions = false;
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
      size_t const numberOfFields( tunnelPath.NumberOfFields() );
      std::vector< double > falseConfiguration( numberOfFields );
      tunnelPath.PutOnPathAt( falseConfiguration,
                              auxiliaryAtRadialInfinity );
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
           fieldIndex < numberOfFields;
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
    }
  }

} /* namespace VevaciousPlusPlus */
