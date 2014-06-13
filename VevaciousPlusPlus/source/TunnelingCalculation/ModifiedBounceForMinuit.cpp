/*
 * ModifiedBounceForMinuit.cpp
 *
 *  Created on: May 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  ModifiedBounceForMinuit::ModifiedBounceForMinuit(
                                    PotentialFunction const& potentialFunction,
                                         size_t const numberOfVaryingPathNodes,
                                       size_t const numberOfSplinesInPotential,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                double const falseVacuumEvaporationTemperature,
                                      size_t const undershootOvershootAttempts,
                                   size_t const maximumMultipleOfLongestLength,
                                  double const initialFractionOfShortestLength,
                                              double const minimumScaleSquared,
                                  double const shootingCloseEnoughThreshold ) :
    ROOT::Minuit2::FCNBase(),
    potentialFunction( potentialFunction ),
    numberOfFields( potentialFunction.NumberOfFieldVariables() ),
    referenceFieldIndex( 0 ),
    pathFromNodes( numberOfFields,
                   referenceFieldIndex,
                   numberOfVaryingPathNodes ),
    numberOfSplinesInPotential( numberOfSplinesInPotential ),
    falseVacuum( falseVacuum ),
    trueVacuum( trueVacuum ),
    zeroTemperatureStraightPath( numberOfFields ),
    zeroTemperatureStraightPathInverseLengthSquared( 0.0 ),
    falseVacuumPotential( potentialFunction( falseVacuum.FieldConfiguration(),
                                             0.0 ) ),
    trueVacuumPotential( potentialFunction( trueVacuum.FieldConfiguration(),
                                            0.0 ) ),
    falseVacuumEvaporationTemperature( falseVacuumEvaporationTemperature ),
    tunnelingScaleSquared( std::max( minimumScaleSquared,
                            potentialFunction.ScaleSquaredRelevantToTunneling(
                                                                   falseVacuum,
                                                              trueVacuum ) ) ),
    shortestLength( 1.0 ),
    longestLength( 1.0 ),
    undershootOvershootAttempts( undershootOvershootAttempts ),
    initialFractionOfShortestLength( initialFractionOfShortestLength ),
    shootingThreshold( shootingCloseEnoughThreshold )
  {
    // We pick the field that has the largest difference between the vacua as
    // the reference field that goes linearly with the auxiliary variable.
    std::vector< double > const&
    falseFieldConfiguration( falseVacuum.FieldConfiguration() );
    std::vector< double > const&
    trueFieldConfiguration( trueVacuum.FieldConfiguration() );
    double greatestFieldDifference( trueFieldConfiguration.front()
                                    - falseFieldConfiguration.front() );
    double currentFieldDifference( greatestFieldDifference );
    zeroTemperatureStraightPath.front() = currentFieldDifference;
    if( greatestFieldDifference < 0.0 )
    {
      greatestFieldDifference = -greatestFieldDifference;
    }
    for( size_t fieldIndex( 1 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentFieldDifference = ( trueFieldConfiguration[ fieldIndex ]
                                 - falseFieldConfiguration[ fieldIndex ] );
      zeroTemperatureStraightPath[ fieldIndex ] = currentFieldDifference;
      zeroTemperatureStraightPathInverseLengthSquared
      += ( currentFieldDifference * currentFieldDifference );
      if( currentFieldDifference < 0.0 )
      {
        currentFieldDifference = -currentFieldDifference;
      }
      if( currentFieldDifference > greatestFieldDifference )
      {
        referenceFieldIndex = fieldIndex;
        greatestFieldDifference = currentFieldDifference;
      }
    }
    zeroTemperatureStraightPathInverseLengthSquared
    = ( 1.0 / zeroTemperatureStraightPathInverseLengthSquared );
    pathFromNodes.SetReferenceField( referenceFieldIndex );
    double lowestScaleSquared( tunnelingScaleSquared );
    double highestScaleSquared( tunnelingScaleSquared );
    double candidateScaleSquared( std::max( minimumScaleSquared,
                                            falseVacuum.LengthSquared() ) );
    if( candidateScaleSquared < lowestScaleSquared )
    {
      lowestScaleSquared = candidateScaleSquared;
    }
    if( candidateScaleSquared > highestScaleSquared )
    {
      highestScaleSquared = candidateScaleSquared;
    }
    candidateScaleSquared = std::max( minimumScaleSquared,
                                      trueVacuum.LengthSquared() );
    if( candidateScaleSquared < lowestScaleSquared )
    {
      lowestScaleSquared = candidateScaleSquared;
    }
    if( candidateScaleSquared > highestScaleSquared )
    {
      highestScaleSquared = candidateScaleSquared;
    }
    shortestLength = ( 1.0 / sqrt( highestScaleSquared ) );
    longestLength = ( 1.0 / sqrt( lowestScaleSquared ) );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "tunnelingScaleSquared = " << tunnelingScaleSquared;
    std::cout << std::endl;/**/
  }

  ModifiedBounceForMinuit::~ModifiedBounceForMinuit()
  {
    // This does nothing.
  }


  // This turns pathParameterization into a PathFieldsAndPotential, by first
  // checking for a non-zero temperature, then setting up the straight path
  // in field space, and projecting the nodes extracted from
  // pathParameterization onto planes perpendicular to the straight path. A
  // few extra bits of information to do with the temperature are also
  // recorded in the PathFieldsAndPotential object returned.
  PathFieldsAndPotential ModifiedBounceForMinuit::DecodePathParameters(
                      std::vector< double > const& pathParameterization ) const
  {
    if( ZeroTemperatureParameterization( pathParameterization ) )
    {
      return pathFromNodes( pathParameterization,
                            zeroTemperatureStraightPath,
                            zeroTemperatureStraightPathInverseLengthSquared,
                            falseVacuum.FieldConfiguration(),
                            falseVacuumPotential,
                            trueVacuumPotential,
                            0.0 );
    }
    double const givenTemperature( pathParameterization.back() );
    PotentialForMinuit potentialForMinuit( potentialFunction );
    potentialForMinuit.SetTemperature( givenTemperature );
    MinuitManager thermalDsbFinder( potentialForMinuit );

    MinuitMinimum
    thermalTrueVacuum( thermalDsbFinder( trueVacuum.FieldConfiguration() ) );
    double const thermalTrueVacuumPotential( thermalTrueVacuum.FunctionValue()
                                     + potentialForMinuit.FunctionAtOrigin() );
    // We undo the offset from potentialFromMinuit.

    double straightPathLengthSquared( 0.0 );
    std::vector< double >&
    straightPath( thermalTrueVacuum.VariableValues() );
    // Might as well put stuff straight into
    // thermalTrueVacuum.VariableValues() as it won't be used later anyway.
    bool const falseMinimumNotEvaporated( givenTemperature
                                         < falseVacuumEvaporationTemperature );
    if( falseMinimumNotEvaporated )
    {
      MinuitMinimum thermalFalseVacuum( thermalDsbFinder(
                                          falseVacuum.FieldConfiguration() ) );
      // We undo the offset from potentialFromMinuit.
      for( unsigned int fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        straightPath[ fieldIndex ]
        -= thermalFalseVacuum.VariableValues()[ fieldIndex ];
        straightPathLengthSquared += ( straightPath[ fieldIndex ]
                                       * straightPath[ fieldIndex ] );
      }
      return pathFromNodes( pathParameterization,
                            straightPath,
                            ( 1.0 / straightPathLengthSquared ),
                            thermalFalseVacuum.VariableValues(),
                            ( thermalFalseVacuum.FunctionValue()
                              + potentialForMinuit.FunctionAtOrigin() ),
                            thermalTrueVacuumPotential,
                            givenTemperature );
    }
    else
    {
      for( unsigned int fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        straightPathLengthSquared += ( straightPath[ fieldIndex ]
                                       * straightPath[ fieldIndex ] );
      }
      return pathFromNodes( pathParameterization,
                            straightPath,
                            ( 1.0 / straightPathLengthSquared ),
                            potentialFunction.FieldValuesOrigin(),
                      potentialFunction( potentialFunction.FieldValuesOrigin(),
                                         givenTemperature ),
                            thermalTrueVacuumPotential,
                            givenTemperature );
    }
  }

  // This returns a polynomial approximation of the potential along the path
  // given by splineCoefficients.
  void ModifiedBounceForMinuit::PotentialAlongPath(
                         PathFieldsAndPotential& pathFieldsAndPotential ) const
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "ModifiedBounceForMinuit::PotentialAlongPath(pathFieldsAndPotential="
    << pathFieldsAndPotential.AsDebuggingString()
    << ") called. numberOfSplinesInPotential = " << numberOfSplinesInPotential;
    std::cout << std::endl;/**/
    std::vector< double > fieldConfiguration( numberOfFields );
    // We choose to take the cos of a linear distribution of equally-spaced
    // points so that we have a finer resolution of the potential near the
    // minima
    double auxiliaryValue( 0.0 );
    for( size_t splinePoint( 1 );
         splinePoint < numberOfSplinesInPotential;
         ++splinePoint )
    {
      auxiliaryValue = ( 0.5 * ( 1.0 - cos( ( (double)splinePoint * M_PI )
                                    / (double)numberOfSplinesInPotential ) ) );
      std::vector< double > const&
      fieldConfiguration( pathFieldsAndPotential.FieldConfiguration(
                                                            auxiliaryValue ) );
      pathFieldsAndPotential.PotentialApproximation().AddPoint( auxiliaryValue,
                                       ( potentialFunction( fieldConfiguration,
                                    pathFieldsAndPotential.GivenTemperature() )
                                         - falseVacuumPotential ) );
    }
    pathFieldsAndPotential.PotentialApproximation().SetSpline(
                                  trueVacuumPotential - falseVacuumPotential );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "ModifiedBounceForMinuit::PotentialAlongPath([pathFieldsAndPotential])"
    << " set pathFieldsAndPotential to be" << std::endl
    << pathFieldsAndPotential.AsDebuggingString();
    std::cout << std::endl;/**/
  }

  // This returns the effective bounce action [S_4 or ((S_3(T)/T + ln(S_3(T)))]
  // calculated under the thin-wall approximation. Before returning, it sets
  // thinWallIsGoodApproximation to be true or false depending on whether
  // the thin-wall approximation is a good approximation or not.
  double ModifiedBounceForMinuit::ThinWallApproximation(
                           PathFieldsAndPotential const pathFieldsAndPotential,
                                      bool& thinWallIsGoodApproximation ) const
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Must remember to test thin-wall approximation!";
    std::cout << std::endl;/**/

    double const potentialDifference( pathFieldsAndPotential.FalsePotential()
                                    - pathFieldsAndPotential.TruePotential() );
    if( potentialDifference
        > ( 1.0E-1 * tunnelingScaleSquared * tunnelingScaleSquared ) )
    {
      thinWallIsGoodApproximation = false;
      return NAN;
    }

    double const oneOverPotentialDifference( 1.0 / potentialDifference );

    // Following Coleman and adapting to thermal tunneling:
    // Quantum:
    // S_E = pi^2 * R^3 * S_1 - 0.5 * pi^2 * R^4 * epsilon
    // R = 3 S_1 / epsilon
    // S_E = -0.5 pi^2 R^3 S_1 = -13.5 pi^2 S_1^4 epsilon^-3
    // Thermal:
    // S_E = 2 * pi * R^2 * S_1 - 4/3 * pi * R^3 * epsilon
    // R = S_1 / epsilon
    // S_E = 2/3 pi R^2 S_1 = 2/3 S_1^3 epsilon^-2

    double integratedAction( 0.0 );
    double potentialValue( 0.0 );
    double const thinWallStep( 0.05 );
    // The start and end of the integral are neglected as they are proportional
    // to the derivatives of the fields with respect to the auxiliary variable,
    // which are zero at the ends of the integral.
    // We also integrate V(field[referenceFieldIndex])^(-1/2) over
    // field[referenceFieldIndex] to estimate the wall thickness, which is then
    // used for comparison to the bubble radius to validate the thin-wall
    // approximation. Because it will be compared to R / 100, accuracy doesn't
    // matter so much, so midpoint summation should suffice.
    double wallThickness( 0.0 );
    for( double thinWallAuxiliary( thinWallStep );
         thinWallAuxiliary < 1.0;
         thinWallAuxiliary += thinWallStep )
    {
      potentialValue
      = pathFieldsAndPotential.PotentialApproximation( thinWallAuxiliary );
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "thinWallAuxiliary = " << thinWallAuxiliary
      << ", potentialValue = " << potentialValue;
      std::cout << std::endl;/**/
      if( potentialValue > 0.0 )
      {
        integratedAction
        += sqrt( pathFieldsAndPotential.FieldDerivativesSquared(
                                        thinWallAuxiliary ) * potentialValue );
        wallThickness += ( thinWallStep / sqrt( potentialValue ) );
      }
    }
    integratedAction *= ( M_SQRT2 * thinWallStep );
    wallThickness *= pathFieldsAndPotential.FieldDerivatives()[
                             referenceFieldIndex ].CoefficientVector().front();
    // This last *= accounts for the change in integration variable to the
    // field linear in a.

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "integratedAction = " << integratedAction
    << ", oneOverPotentialDifference = " << oneOverPotentialDifference
    << ", wallThickness = " << wallThickness;
    std::cout << std::endl;/**/


    // We return the thin-wall approximation of the bounce action only if the
    // bubble radius is sufficiently large compared to the tunneling scale,
    // which is a good indicator for the validity of the approximation.
    if( ( integratedAction * oneOverPotentialDifference )
        < ( 1.0E+2 * wallThickness ) )
    {
      thinWallIsGoodApproximation = false;
      return NAN;
    }
    thinWallIsGoodApproximation = false;
    if( pathFieldsAndPotential.NonZeroTemperature() )
    {
      double const bounceAction( ( 2.0 * M_PI
                       * integratedAction * integratedAction * integratedAction
                    * oneOverPotentialDifference * oneOverPotentialDifference )
                    / pathFieldsAndPotential.GivenTemperature() );
      // This is at least 10^4 * 2 pi * ( S_1 / ( Q^2 T ) ), so S_3(T)/T is
      // likely to be at least 10^3, leading to a totally negligible
      // tunneling probability...
      return ( ( bounceAction / pathFieldsAndPotential.GivenTemperature() )
               + log( bounceAction ) );
    }
    else
    {
      return ( -13.5 * M_PI * M_PI
                     * integratedAction * integratedAction
                     * integratedAction * integratedAction
                     * oneOverPotentialDifference
                     * oneOverPotentialDifference
                     * oneOverPotentialDifference );
      // This is at least 3^4 * 10^4 * 13.5 pi^2 * ( S_1 / ( Q^3 ) ), so
      // S_4 is likely to be at least 10^6, leading to a totally negligible
      // tunneling probability...
    }
  }

  // This sets up the bubble profile, numerically integrates the bounce action
  // over it, and then returns effective bounce action
  // [S_4 or ((S_3(T)/T + ln(S_3(T)))].
  double ModifiedBounceForMinuit::EffectiveBounceAction(
                   PathFieldsAndPotential const& pathFieldsAndPotential ) const
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "ModifiedBounceForMinuit::EffectiveBounceAction("
    << " pathFieldsAndPotential =" << std::endl
    << pathFieldsAndPotential.AsDebuggingString() << " ) called.";
    std::cout << std::endl;/**/
    BubbleProfile bubbleProfile( pathFieldsAndPotential,
                          ( initialFractionOfShortestLength * shortestLength ),
                                 longestLength );
    std::vector< BubbleRadialValueDescription > const&
    auxiliaryProfile( bubbleProfile.DampedProfile( undershootOvershootAttempts,
                                                   shootingThreshold ) );
    // The bounce action density at r = 0 is 0 by merit of
    // r_0^dampingFactor = 0,
    // dp/dr = 0 at r = 0,
    // and potentialApproximation(p(0)) = 0 by construction.
    // The smallest r recorded by bubbleProfile is non-zero because of this.
    // Thus we add a sphere of 3 or 4 dimensions with just the potential
    // contribution, assuming the kinetic contribution up to the 1st radius is
    // negligible, before adding in spherical layers of integrand. (The solid
    // angle factor is left until afterwards.)
    double currentAuxiliary( auxiliaryProfile.front().auxiliaryValue );
    double currentRadius( auxiliaryProfile.front().radialValue );
    double currentIntegrand( pathFieldsAndPotential.PotentialApproximation(
                                                             currentAuxiliary )
                             * currentRadius * currentRadius );
    double previousIntegrand( NAN );
    double bounceAction( NAN );
    if( pathFieldsAndPotential.NonZeroTemperature() )
    {
      bounceAction = ( currentIntegrand / 3.0 );
    }
    else
    {
      currentIntegrand *= currentRadius;
      bounceAction = ( 0.25 * currentIntegrand );
    }
    double previousRadius( currentIntegrand );
    double kineticTerm( NAN );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "keeping track of parts of integral. pathFieldsAndPotential ="
    << std::endl << pathFieldsAndPotential.AsDebuggingString();
    double debugKineticSum( 0.0 );
    double debugPotentialSum( 0.0 );
    std::cout << std::endl;/**/
    for( size_t radiusIndex( 1 );
         radiusIndex < auxiliaryProfile.size();
         ++radiusIndex )
    {
      previousRadius = currentRadius;
      previousIntegrand = currentIntegrand;
      currentRadius = auxiliaryProfile[ radiusIndex ].radialValue;
      currentAuxiliary = auxiliaryProfile[ radiusIndex ].auxiliaryValue;
      kineticTerm = auxiliaryProfile[ radiusIndex ].auxiliarySlope;
      kineticTerm
      *= ( 0.5 * kineticTerm
        * pathFieldsAndPotential.FieldDerivativesSquared( currentAuxiliary ) );
      currentIntegrand
      = ( ( kineticTerm
          + pathFieldsAndPotential.PotentialApproximation( currentAuxiliary ) )
                           * currentRadius * currentRadius );
      if( pathFieldsAndPotential.NonZeroTemperature() )
      {
        currentIntegrand *= currentRadius;
      }
      bounceAction += ( ( currentIntegrand + previousIntegrand )
                        * ( currentRadius - previousRadius ) );
      // A common factor of 1/2 is being left until after the loop.
      // debugging:
      /**/
      debugKineticSum
      += ( kineticTerm
          * currentRadius * currentRadius * currentRadius
          * ( currentRadius - previousRadius ) );
      debugPotentialSum
      += ( pathFieldsAndPotential.PotentialApproximation( currentAuxiliary )
           * currentRadius * currentRadius * currentRadius
           * ( currentRadius - previousRadius ) );
      std::cout << std::endl << "debugging:"
      << std::endl
      << "radiusIndex = " << radiusIndex
      << ", r = " << currentRadius
      << ", dr = " << ( currentRadius - previousRadius )
      << ", p = " << currentAuxiliary
      << ", |dp/dr|^2 = "
      << pow( auxiliaryProfile[ radiusIndex ].auxiliarySlope,
              2 )
      << ", df/dp.df/dp = "
      << pathFieldsAndPotential.FieldDerivativesSquared( currentAuxiliary )
      << ", V(p) = "
      << pathFieldsAndPotential.PotentialApproximation( currentAuxiliary )
      << ", debugKineticSum = " << debugKineticSum
      << ", debugPotentialSum = " << debugPotentialSum;
      std::cout << std::endl;/**/
    }
    // The common factor of 1/2 is combined with the solid angle of
    // 2 pi^2 (quantum) or 4 pi (thermal):
    if( pathFieldsAndPotential.NonZeroTemperature() )
    {
      bounceAction *= ( 2.0 * M_PI );
    }
    else
    {
      bounceAction *= ( M_PI * M_PI );
    }

    if( pathFieldsAndPotential.NonZeroTemperature() )
    {
      return ( ( bounceAction / pathFieldsAndPotential.GivenTemperature() )
               + log( bounceAction ) );
    }
    else
    {
      return bounceAction;
    }
  }

} /* namespace VevaciousPlusPlus */
