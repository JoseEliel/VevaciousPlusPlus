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
                                      size_t const potentialApproximationPower,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                double const falseVacuumEvaporationTemperature,
                                      size_t const undershootOvershootAttempts,
                                   size_t const maximumMultipleOfLongestLength,
                                  double const initialFractionOfShortestLength,
                               size_t const energyConservingUndershootAttempts,
                                              double const minimumScaleSquared,
                                  double const shootingCloseEnoughThreshold ) :
    ROOT::Minuit2::FCNBase(),
    potentialFunction( potentialFunction ),
    numberOfFields( potentialFunction.NumberOfFieldVariables() ),
    referenceFieldIndex( 0 ),
    pathFromNodes( numberOfFields,
                   referenceFieldIndex,
                   numberOfVaryingPathNodes ),
    potentialApproximationPower( potentialApproximationPower ),
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
    energyConservingUndershootAttempts( energyConservingUndershootAttempts ),
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
    std::vector< double > fieldConfiguration( numberOfFields );
    // Here we set up the linear system to solve for the coefficients of the
    // polynomial approximation of the potential:
    size_t const
    numberOfNonZeroCoefficients( potentialApproximationPower - 2 );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "ModifiedBounceForMinuit::PotentialAlongPath(...) called."
    << " numberOfNonZeroCoefficients = " << numberOfNonZeroCoefficients;
    std::cout << std::endl;/**/
    // The potential is approximated as a polynomial with zero coefficients
    // (for a^0 and a^1, so that a = 0 is an extremum) which will be ignored
    // for inverting the matrix equation. Furthermore, the final coefficient is
    // related to the other coefficients by the condition that the derivative
    // is 0 for a = 1, so for n = potentialApproximationPower
    // -n c_n = 2 c_2 + 3 c_3 + ...
    // making the equation defining the potential approximation become
    // V(a) = c_2 ( a^2 - (2/n) a^n ) + c_3 ( a^3 - (3/n) a^n ) + ...
    // and c_n will be explicitly put into the approximation object afterwards.
    double const
    auxiliaryStep( 1.0 / (double)numberOfNonZeroCoefficients );
    double auxiliaryValue( 0.0 );
    Eigen::VectorXd potentialValues( numberOfNonZeroCoefficients );
    Eigen::MatrixXd coefficientMatrix( numberOfNonZeroCoefficients,
                                       numberOfNonZeroCoefficients );
    // We want to solve coefficientMatrix * c = potentialValues for the
    // coefficients c. Since the potential approximation begins with the
    // quadratic term, coefficientMatrix, for n auxiliary values x_1 to x_n and
    // potentialApproximationPower = d, looks like
    // [ [ x_1^2, x_1^3, ..., x_1^d ],
    //   [ x_2^2, x_2^3, ..., x_2^d ],
    //   ...,
    //   [ x_n^2, x_n^3, ..., x_n^d ] ], and is n * (d-1) in size. However, we
    // also fix c_d = (-1/d) * (2 c_2 + 3 c_3 + ... + (d-1) c_(d-1)) so that
    // V(1) is also a minimum (well, extremum, but we hope that the
    // approximation works out so that it is a minimum). Hence we can just take
    // coefficientMatrix to be n * (d-2) in size, noting that c_j multiplies
    // ( x^j - (j/d) x^d ) having substituted in for c_d. Therefore,
    // coefficientMatrix looks like
    // [ [ ( x_1^2 - (2/d) x_1^d ), ( x_1^3- (3/d) x_1^d ), ...,
    //                                        ( x_1^(d-1)- ((d-1)/d) x_1^d ) ],
    //   ...,
    //   [ ( x_n^2 - (2/d) x_n^d ), ( x_n^3- (3/d) x_n^d ), ...,
    //                                       ( x_n^(d-1)- ((d-1)/d) x_n^d ) ] ]
    // and also x_j is simply (j/n).

    double finalCoefficientPart( NAN );
    double const finalCoefficientDenominatorInverse( 1.0 /
                                            ( pow( numberOfNonZeroCoefficients,
                                                  potentialApproximationPower )
                                             * potentialApproximationPower ) );
    // This gives the recurring (1/d(n^d)) from (j/d) x_k^d = (j/d) (k/n)^d.
    // It is multiplied by j k^d.
    // This corresponds now to k = whichStep + 1, j = whichPower in the
    // following nested loops.
    for( unsigned int whichStep( 0 );
         whichStep < ( numberOfNonZeroCoefficients - 1 );
         ++whichStep )
    {
      auxiliaryValue += auxiliaryStep;
      // Now auxiliaryValue = x_k = k/n.
      std::vector< double > const&
      fieldConfiguration( pathFieldsAndPotential.FieldConfiguration(
                                                            auxiliaryValue ) );
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "fieldConfiguration = { ";
      for( std::vector< double >::const_iterator
           fieldValue( fieldConfiguration.begin() );
           fieldValue < fieldConfiguration.end();
           ++fieldValue )
      {
        if( fieldValue != fieldConfiguration.begin() )
        {
          std::cout << ", ";
        }
        std::cout << *fieldValue;
      }
      std::cout << " }. Going to set potentialValues( " << whichStep
      << " ) to be " << ( potentialFunction( fieldConfiguration,
                                    pathFieldsAndPotential.GivenTemperature() )
                          - falseVacuumPotential );
      std::cout << std::endl;/**/
      potentialValues( whichStep ) = ( potentialFunction( fieldConfiguration,
                                    pathFieldsAndPotential.GivenTemperature() )
                                       - falseVacuumPotential );

      finalCoefficientPart = ( pow( ( whichStep + 1 ),
                                    potentialApproximationPower )
                               * finalCoefficientDenominatorInverse );
      // Now finalCoefficientPart = (k^d)/(d(n^d)) = ( (j/d) x_k^d )/j.
      for( unsigned int whichPower( 2 );
           whichPower < potentialApproximationPower;
           ++whichPower )
      {
        coefficientMatrix( whichStep,
                           ( whichPower - 2 ) )
        = ( pow( auxiliaryValue,
                 whichPower )
            - ( (double)whichPower * finalCoefficientPart ) );
        // Now coefficientMatrix( k, (j-2) ) = x_k^j - (j/d) x_k^d.
      }
    }
    potentialValues( numberOfNonZeroCoefficients - 1 )
    = ( trueVacuumPotential - falseVacuumPotential );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl;
    std::cout << "potentialValues =" << std::endl << potentialValues;
    std::cout << std::endl;/**/
    finalCoefficientPart = ( 1.0 / potentialApproximationPower );
    for( unsigned int whichPower( 0 );
         whichPower < numberOfNonZeroCoefficients;
         ++whichPower )
    {
      coefficientMatrix( ( numberOfNonZeroCoefficients - 1 ),
                         whichPower )
      = ( 1.0 - ( (double)( whichPower + 2 ) * finalCoefficientPart ) );
      // Now coefficientMatrix( n, (j-2) ) = 1^j - (j/d) 1^d.
    }
    std::cout << std::endl << "coefficientMatrix =" << std::endl
    << coefficientMatrix;

    pathFieldsAndPotential.SetPotential(
            coefficientMatrix.colPivHouseholderQr().solve( potentialValues ) );
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
        > ( 1.0E-3 * tunnelingScaleSquared * tunnelingScaleSquared ) )
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
    double fieldDerivative( 0.0 );
    double derivativeSquared( 0.0 );
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
      derivativeSquared = 0.0;
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        fieldDerivative = pathFieldsAndPotential.FieldDerivative( fieldIndex,
                                                           thinWallAuxiliary );
        derivativeSquared += ( fieldDerivative * fieldDerivative );
      }
      integratedAction += sqrt( derivativeSquared * potentialValue );
      wallThickness += ( thinWallStep / sqrt( potentialValue ) );
    }
    integratedAction *= ( M_SQRT2 * thinWallStep );
    wallThickness *= pathFieldsAndPotential.FieldDerivatives()[
                             referenceFieldIndex ].CoefficientVector().front();
    // This last *= accounts for the change in integration variable to the
    // field linear in a.


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
    bubbleProfile.UndampedUndershoot( energyConservingUndershootAttempts );
    std::vector< BubbleRadialValueDescription > const&
    auxiliaryProfile( bubbleProfile.DampedProfile( undershootOvershootAttempts,
                                                   shootingThreshold ) );
    double bounceAction( 0.0 );
    double previousIntegrand( 0.0 );
    // The bounce action density at r = 0 is 0 by merit of
    // r_0^dampingFactor = 0,
    // dp/dr = 0 at r = 0,
    // and potentialApproximation(p(0)) = 0 by construction.
    // The smallest r recorded by bubbleProfile is non-zero because of this.
    double currentAuxiliary( NAN );
    double currentIntegrand( 0.0 );
    double previousRadius( 0.0 );
    double currentRadius( 0.0 );
    double kineticTerm( NAN );
    for( unsigned int radiusIndex( 0 );
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
