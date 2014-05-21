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
                                unsigned int const potentialApproximationPower,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                             double const falseVacuumEvaporationTemperature ) :
    ROOT::Minuit2::FCNBase(),
    potentialFunction( potentialFunction ),
    numberOfFields( potentialFunction.NumberOfFieldVariables() ),
    potentialApproximationPower( potentialApproximationPower ),
    falseVacuum( falseVacuum ),
    trueVacuum( trueVacuum ),
    falseVacuumEvaporationTemperature( falseVacuumEvaporationTemperature )
  {
    // This constructor is just an initialization list.
  }

  ModifiedBounceForMinuit::~ModifiedBounceForMinuit()
  {
    // This does nothing.
  }


  // This implements operator() for FCNBase, the function that MINUIT will
  // minimize. The values of splineCoefficients should be sets of n
  // coefficients for polynomials for each of the n fields, plus the final
  // element of splineCoefficients should be the temperature.
  double ModifiedBounceForMinuit::operator()(
                        std::vector< double > const& splineCoefficients ) const
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "ModifiedBounceForMinuit::operator() doesn't work at all at the moment."
    << " We need to think a lot harder about how to implement it.";
    std::cout << std::endl;
    return 0.0;/**/


    // Current outline (code below doesn't do most of this yet):
    // 1) Take splineCoefficients as proposed dependence of field configuration
    //    f on auxiliary variable a. An extra power of a is given with a
    //    coefficient to bring f to the true vacuum at a = 1, as
    //    splineCoefficients alone does not have a fixed endpoint at a = 1 for
    //    different splineCoefficients input.
    // 2) Evaluate the potential at ( potentialApproximationPower - 2 ) values
    //    for a between a = 0 and a = 1 inclusive, giving V(a) as a polynomial
    //    in a (though V(a) is an approximation, f(a) is exact by
    //    construction). (The number of points to evaluate V(a) would be
    //    potentialApproximationPower + 1 for potentialApproximationPower
    //    coefficients of non-zero powers of a plus a constant coefficient, but
    //    2 of the coefficients are fixed by the requirement that dV/da is 0 at
    //    a = 0 and a = 1, and the constant coefficient by taking V(0) to be 0,
    //    which also is convenient for evaluating the action. Hence V(a) is
    //    actually a polynomial approximation of the potential difference from
    //    the false vacuum along the path given by f(a).)
    // 3a) Return early with a thin-wall result if appropriate.
    // 3b) Find the 1-dimensional bubble profile by solving the radial bubble
    //     equations of motion along the path from a = 0 to a = a_crit, where
    //     a_crit is the value of a for which the bubble is critical, and is
    //     found by the undershoot/overshoot method (0.0 < a_crit < 1.0). This
    //     probably will involve using the Boost library odeint.
    // 4) Combine f(a) (exact) and V(a) (approximate) with a(r) to numerically
    //    integrate the bounce action.
    // Note that this is only really the bounce action (within the validity of
    // the approximations of truncated exapnsions) if the equations of motion
    // perpendicular to f(a) are also satisfied. However, they would be
    // satisfied for a modified potential which increases the energy barrier on
    // the side of the path until the "force" on the path balances, and such a
    // modified potential still has the same true and false vacua, and its
    // decay width is given by this calculated bounce action. The unmodified
    // potential with the possibility for a path with a lower energy barrier
    // must necessarily have at least as large a decay width as the modified
    // potential, hence the calculated bounce action here is an upper bound on
    // the true bounce action, and minimizing this should lead one to the true
    // minimal bounce action.


    double const givenTemperature( splineCoefficients.back() );
    bool const nonZeroTemperature( givenTemperature > 0.0 );

    std::vector< double > falseConfiguration( numberOfFields,
                                              0.0 );
    double falsePotentialValue( 0.0 );
    if( givenTemperature < falseVacuumEvaporationTemperature )
    {
      PotentialForMinuit potentialForMinuit( potentialFunction );
      MinuitManager thermalDsbFinder( potentialForMinuit );
      MinuitMinimum thermalFalseVacuum
      = thermalDsbFinder( falseVacuum.FieldConfiguration() );
      falsePotentialValue = thermalFalseVacuum.FunctionValue();
      falseConfiguration = thermalFalseVacuum.VariableValues();
    }
    else
    {
      falsePotentialValue = potentialFunction( falseConfiguration,
                                               givenTemperature );
    }
    std::vector< double > fieldConfiguration( falseConfiguration );

    // Here we set up the linear system to solve for the coefficients of the
    // polynomial approximation of the potential:
    unsigned int const
    numberOfNonZeroCoefficients( potentialApproximationPower - 2 );
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
    double const
    finalCoefficientPart( pow( auxiliaryValue,
                               potentialApproximationPower )
                          / (double)potentialApproximationPower );
    for( unsigned int whichStep( 0 );
         whichStep < ( numberOfNonZeroCoefficients - 1 );
         ++whichStep )
    {
      auxiliaryValue += auxiliaryStep;
      ConfigurationFromSplines( fieldConfiguration,
                                splineCoefficients,
                                auxiliaryValue );
      potentialValues( whichStep ) = ( potentialFunction( fieldConfiguration )
                                       - falsePotentialValue );

      for( unsigned int whichPower( 2 );
           whichPower < potentialApproximationPower;
           ++whichPower )
      {
        coefficientMatrix( whichStep,
                           ( whichPower - 2 ) )
        = ( pow( auxiliaryValue,
                 whichPower )
            - ( (double)whichPower * finalCoefficientPart ) );
      }
    }
    potentialValues( numberOfNonZeroCoefficients - 1 )
    = trueVacuum.PotentialValue();
    for( unsigned int whichPower( 0 );
         whichPower < numberOfNonZeroCoefficients;
         ++whichPower )
    {
      coefficientMatrix( ( numberOfNonZeroCoefficients - 1 ),
                         whichPower ) = 1.0;
    }

    SimplePolynomial
    potentialApproximation( coefficientMatrix.colPivHouseholderQr().solve(
                                                             potentialValues ),
                            2 );
    std::vector< double >&
    potentialVector( potentialApproximation.CoefficientVector() );
    double finalCoefficientTimesMaxPower( 0.0 );
    for( unsigned int whichPower( 2 );
         whichPower < potentialApproximationPower;
         ++whichPower )
    {
      finalCoefficientTimesMaxPower
      -= ( (double)whichPower * potentialVector[ whichPower ] );
    }
    potentialVector.push_back( finalCoefficientTimesMaxPower
                               / (double)potentialApproximationPower );
    double const
    tunnelingScaleSquared( potentialFunction.ScaleSquaredRelevantToTunneling(
                                                                   falseVacuum,
                                                                trueVacuum ) );
    double const tunnelingScale( sqrt( tunnelingScaleSquared ) );
    std::vector< SimplePolynomial > fieldDerivatives;
    FieldDerivatives( splineCoefficients,
                      fieldDerivatives );
    // We return a thin-wall approximation if appropriate:
    if( ( ( falseVacuum.PotentialValue() - trueVacuum.PotentialValue() )
          / ( tunnelingScaleSquared * tunnelingScaleSquared ) ) < 1.0E-3 )
    {
      double const
      oneOverPotentialDifference( 1.0 / ( falseVacuum.PotentialValue()
                                          - trueVacuum.PotentialValue() ) );

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
      double fieldDerivative( 0.0 );
      double derivativeSquared( 0.0 );
      double const thinWallStep( 0.05 );
      // The start and end of the integral are neglected as they are
      // proportional to the derivatives of the fields with respect to the
      // auxiliary variable, which are zero at the ends of the integral.
      for( double thinWallAuxiliary( thinWallStep );
           thinWallAuxiliary < 1.0;
           thinWallAuxiliary += thinWallStep )
      {
        derivativeSquared = 0.0;
        for( size_t fieldIndex( 0 );
             fieldIndex < numberOfFields;
             ++fieldIndex )
        {
          fieldDerivative = fieldDerivatives( thinWallAuxiliary );
          derivativeSquared += ( fieldDerivative * fieldDerivative );
        }
        integratedAction += sqrt( derivativeSquared
                               * potentialApproximation( thinWallAuxiliary ) );
      }
      integratedAction *= ( M_SQRT2 * thinWallStep );

      // We return the thin-wall approximation of the bounce action only if the
      // bubble radius is sufficiently large compared to the tunneling scale,
      // which is a good indicator for the validity of the approximation.
      // (The comparison here is that the thermal bubble radius > 100 Q, or the
      // T=0 radius > 300 Q.)
      if( ( integratedAction * oneOverPotentialDifference * tunnelingScale )
          > 1.0E+2 )
      {
        if( nonZeroTemperature )
        {
          return ( ( 2.0 * M_PI * integratedAction * integratedAction
                     * integratedAction * oneOverPotentialDifference
                     * oneOverPotentialDifference ) / givenTemperature );
          // This is at least 10^4 * 2 pi * ( S_1 / ( Q^2 T ) ), so S_3(T)/T is
          // likely to be at least 10^3, leading to a totally negligible
          // tunneling probability...
        }
        else
        {
          return ( -13.5 * M_PI * M_PI * integratedAction * integratedAction
                   * integratedAction * integratedAction
                   * oneOverPotentialDifference * oneOverPotentialDifference
                   * oneOverPotentialDifference );
          // This is at least 3^4 * 10^4 * 13.5 pi^2 * ( S_1 / ( Q^3 ) ), so
          // S_4 is likely to be at least 10^6, leading to a totally negligible
          // tunneling probability...
        }
      }
      // If the thin-wall radius wasn't large enough, nothing is returned yet,
      // and we go on to try undershooting/overshooting.
    }
    unsigned int dampingFactor( 4 );
    if( nonZeroTemperature )
    {
      dampingFactor = 3;
    }
    BubbleProfiler bubbleProfiler( potentialApproximation,
                                   fieldDerivatives,
                                   dampingFactor );
    OdeintBubbleObserver odeintBubbleObserver;

    // Now we have to find the value of the auxiliary variable which is closest
    // to the perfect value that neither undershoots nor overshoots as the
    // radial variable goes to infinity. Unfortunately we have to cope with
    // finite calculations though...
    // The strategy is to start with definite undershoot a (largest a such that
    // V(a) > 0.0, decided from potentialValues) and a definite overshoot
    // (a = 1.0), and narrow the range with definite undershoots and overshoots
    // for a given maximum radial value. The maximum radial value is increased
    // and an undecided a is found each time until the a at maximum radial
    // value is < 1.0E-3 or something.
    // HERE! REMEMBER TO RESET odeintBubbleObserver IN LOOP!



    double solidAngle( 4.0 * M_PI );
    if( nonZeroTemperature )
    {
      solidAngle = ( 2.0 * M_PI * M_PI );
    }

    double fieldGradient( 0.0 );
    double gradientDotGradient( 0.0 );

    double radiusValue( 0.0 );
    double minusRadiusDerivative( 0.0 );
    double factorTimesVolume( 0.0 );

    std::vector< double > falseConfiguration( numberOfFields,
                                              0.0 );
    if( givenTemperature < falseVacuumEvaporationTemperature )
    {
      PotentialForMinuit potentialForMinuit( potentialFunction );
      MinuitManager thermalDsbFinder( potentialForMinuit );
      falseConfiguration
      = thermalDsbFinder( falseVacuum.FieldConfiguration() ).VariableValues();
    }
    double const falsePotential( potentialFunction( falseConfiguration,
                                                    givenTemperature ) );
    std::vector< double > previousFieldConfiguration( falseConfiguration );
    std::vector< double >
    currentFieldConfiguration( previousFieldConfiguration );
    std::vector< double > nextFieldConfiguration( falseConfiguration );
    ConfigurationFromSplines( currentFieldConfiguration,
                              splineCoefficients,
                              auxiliaryValue );

    // The action is evaluated with the composite Simpson's rule, using an
    // auxiliary variable that goes from 0 at the false vacuum at a radius of
    // infinity to 1 at the vacuum at the center of the bubble.

    // The contribution to the bounce action at bubble radius of infinity
    // ( corresponding to the auxiliary variable being zero) should be zero.
    double bounceAction( 0.0 );
    for( unsigned int whichStep( 1 );
         whichStep < numberOfPathIntervals;
         ++whichStep )
    {
      previousFieldConfiguration = currentFieldConfiguration;
      currentFieldConfiguration = nextFieldConfiguration;
      radiusValue = bubbleRadiusFromAuxiliary.RadiusValue( auxiliaryValue );
      minusRadiusDerivative
      = -bubbleRadiusFromAuxiliary.RadiusDerivative( auxiliaryValue );
      // Since the radius should go from infinity to 0 as auxiliaryValue goes
      // from 0 to 1, the derivative is negative.
      factorTimesVolume = ( (double)( 1 + ( whichStep % 2 ) )
                            * solidAngle * radiusValue * radiusValue
                            * minusRadiusDerivative );
      if( !nonZeroTemperature )
      {
        factorTimesVolume *= radiusValue;
      }

      auxiliaryValue += auxiliaryStep;
      // We advance the field configuration to the end of the next interval so
      // that we can take the average derivative of the fields.
      nextFieldConfiguration = falseConfiguration;
      ConfigurationFromSplines( nextFieldConfiguration,
                                splineCoefficients,
                                auxiliaryValue );
      gradientDotGradient = 0.0;
      for( unsigned int fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        fieldGradient = ( ( nextFieldConfiguration[ fieldIndex ]
                            - previousFieldConfiguration[ fieldIndex ] )
                          / ( 2.0 * auxiliaryStep ) );
        gradientDotGradient += ( fieldGradient * fieldGradient );
      }

      bounceAction
      += ( factorTimesVolume
           * ( ( 0.5 * gradientDotGradient )
               + potentialFunction( currentFieldConfiguration,
                                    givenTemperature )
               - falsePotential ) );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "step " << whichStep << ", currentFieldConfiguration = {";
      for( std::vector< double >::const_iterator
           fieldValue( currentFieldConfiguration.begin() );
           fieldValue < currentFieldConfiguration.end();
           ++fieldValue )
      {
        std::cout << " " << *fieldValue;
      }
      std::cout << " }, 0.5 * gradientDotGradient = "
      << ( 0.5 * gradientDotGradient ) << ", potential difference = "
      << ( potentialFunction( currentFieldConfiguration,
                              givenTemperature ) - falsePotential )
      << ", factorTimesVolume = " << factorTimesVolume
      << ", bounceAction = " << bounceAction;
      std::cout << std::endl;/**/
    }
    // There was a common factor of 2.0 taken out of the above, so we account
    // for it now:
    bounceAction += bounceAction;

    // The last value is special, as it corresponds to the middle of the
    // bubble: we force the field derivative to zero.
    factorTimesVolume = ( 0.5 * solidAngle * radiusValue * radiusValue );
    if( !nonZeroTemperature )
    {
      factorTimesVolume *= radiusValue;
    }
    bounceAction += ( factorTimesVolume
                      * ( potentialFunction( currentFieldConfiguration,
                                             givenTemperature )
                          - falsePotential ) );

    // Now we multiply by the common factor:
    bounceAction *= ( auxiliaryStep / 3.0 );

    if( nonZeroTemperature )
    {
      return ( ( bounceAction / givenTemperature ) + log( bounceAction ) );
    }
    else
    {
      return bounceAction;
    }
  }

  void ModifiedBounceForMinuit::FieldsAsSimplePolynomials(
                               std::vector< double > const& splineCoefficients,
                   std::vector< SimplePolynomial >& fieldsAsPolynomials ) const
  {
    unsigned int coefficientsPerField( ( ( splineCoefficients.size() - 2 )
                                       / numberOfFields ) + 2 );
    fieldsAsPolynomials.resize( numberOfFields,
                                SimplePolynomial( coefficientsPerField ) );
    // The last element of splineCoefficients is a temperature, so should not
    // be used. The ( ( size - 2 / numberOfFields ) + 2 ) integer is to ensure
    // that each SimplePolynomial has enough elements. The last coefficient of
    // each field is implicit, and is fixed by the requirement that the field
    // path goes to trueVacuum for a=1.
    for( unsigned int splineIndex( 0 );
         splineIndex < ( splineCoefficients.size() - 1 );
         ++splineIndex )
    {
      fieldsAsPolynomials[ splineIndex
                           % numberOfFields ].CoefficientVector()[ splineIndex
                                                             / numberOfFields ]
      = splineCoefficients[ splineIndex ];
    }
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double coefficientSum( 0.0 );
      for( std::vector< double >::const_iterator coefficientValue(
               fieldsAsPolynomials[ fieldIndex ].CoefficientVector().begin() );
           coefficientValue
           < ( fieldsAsPolynomials[ fieldIndex ].CoefficientVector().end()
               - 1 );
           ++coefficientValue )
      {
        coefficientSum += *coefficientValue;
      }
      fieldsAsPolynomials[ fieldIndex ].CoefficientVector().back()
      = ( trueVacuum.FieldConfiguration()[ fieldIndex ]
          - falseVacuum.FieldConfiguration()[ fieldIndex ]
          - coefficientSum );
      // This ensures that the field configuration for auxiliary variable = 1.0
      // goes to that of the true vacuum.
    }
  }

  // This puts the derivatives of the fields directly into fieldDerivatives
  // from splineCoefficients.
  void ModifiedBounceForMinuit::FieldDerivatives(
                               std::vector< double > const& splineCoefficients,
                      std::vector< SimplePolynomial >& fieldDerivatives ) const
  {
    unsigned int coefficientsPerField( ( ( splineCoefficients.size() - 2 )
                                       / numberOfFields ) + 1 );
    fieldDerivatives.resize( numberOfFields,
                             SimplePolynomial( coefficientsPerField ) );
    // The last element of splineCoefficients is a temperature, so should not
    // be used. The ( ( size - 2 / numberOfFields ) + 1 ) integer is to ensure
    // that each SimplePolynomial has enough elements. The last coefficient of
    // each field is implicit, and is fixed by the requirement that the field
    // path goes to trueVacuum for a=1.
    for( unsigned int splineIndex( 0 );
         splineIndex < ( splineCoefficients.size() - 1 - numberOfFields );
         ++splineIndex )
    {
      fieldDerivatives[ splineIndex
                        % numberOfFields ].CoefficientVector()[ ( splineIndex
                                                           / numberOfFields ) ]
      = splineCoefficients[ splineIndex + numberOfFields ];
    }
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double coefficientSum( 0.0 );
      unsigned int whichPower( 0 );
      while( whichPower < ( coefficientsPerField - 1 ) )
      {
        coefficientSum
        += fieldDerivatives[ fieldIndex ].CoefficientVector()[ whichPower ];
        fieldDerivatives[ fieldIndex ].CoefficientVector()[ whichPower ]
        *= (double)( ++whichPower );
      }
      fieldDerivatives[ fieldIndex ].CoefficientVector().back()
      = ( trueVacuum.FieldConfiguration()[ fieldIndex ]
          - falseVacuum.FieldConfiguration()[ fieldIndex ]
          - coefficientSum
          - splineCoefficients[ fieldIndex ] );
      // This ensures that the field configuration for auxiliary variable = 1.0
      // goes to that of the true vacuum.
      fieldDerivatives[ fieldIndex ].CoefficientVector().back()
       *= (double)coefficientsPerField;
    }
  }

} /* namespace VevaciousPlusPlus */
