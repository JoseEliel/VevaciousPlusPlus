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
                             double const falseVacuumEvaporationTemperature ) :
    ROOT::Minuit2::FCNBase(),
    potentialFunction( potentialFunction ),
    numberOfFields( potentialFunction.NumberOfFieldVariables() ),
    potentialApproximationPower( potentialApproximationPower ),
    falseVacuum( falseVacuum ),
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
    // 2) Evaluate the potential at ( potentialApproximationPower + 1 ) values
    //    for a between a = 0 and a = 1 inclusive, giving V(a) as a polynomial
    //    in a (though V(a) is an approximation, f(a) is exact by
    //    construction).
    // 3) Find the 1-dimensional bubble profile by solving the radial bubble
    //    equations of motion along the path from a = 0 to a = a_crit, where
    //    a_crit is the value of a for which the bubble is critical, and is
    //    found by the undershoot/overshoot method (0.0 < a_crit < 1.0). This
    //    probably will involve using the Boost library odeint.
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
    if( givenTemperature < falseVacuumEvaporationTemperature )
    {
      PotentialForMinuit potentialForMinuit( potentialFunction );
      MinuitManager thermalDsbFinder( potentialForMinuit );
      falseConfiguration
      = thermalDsbFinder( falseVacuum.FieldConfiguration() ).VariableValues();
    }
    double const falsePotential( potentialFunction( falseConfiguration,
                                                    givenTemperature ) );
    std::vector< double > fieldConfiguration( falseConfiguration );

    // Here we set up the linear system to solve for the coefficients of the
    // polynomial approximation of the potential:
    double const auxiliaryStep( 1.0 / (double)potentialApproximationPower );
    double auxiliaryValue( 0.0 );
    Eigen::VectorXd potentialValues( potentialApproximationPower );
    Eigen::MatrixXd coefficientMatrix( potentialApproximationPower,
                                       potentialApproximationPower );
    for( unsigned int whichStep( 0 );
         whichStep < potentialApproximationPower;
         ++whichStep )
    {
      auxiliaryValue += auxiliaryStep;
      ConfigurationFromSplines( fieldConfiguration,
                                splineCoefficients,
                                auxiliaryValue );
      potentialValues( whichStep ) = potentialFunction( fieldConfiguration );

      for( unsigned int whichPower( 0 );
           whichPower < potentialApproximationPower;
           ++whichPower )
      {
        coefficientMatrix( whichStep, whichPower ) = pow( auxiliaryValue,
                                                          ( whichPower + 1 ) );
      }
    }

    SimplePolynomial potentialApproximationOverAuxiliary(
                                 coefficientMatrix.colPivHouseholderQr().solve(
                                                           potentialValues ) );
    // Because the potential is approximated as a polynomial with zero
    // coefficient which was ignored for inverting the matrix equation,
    // potentialApproximationOverAuxiliary is the potential divided by the
    // auxiliary variable, which is convenient for creating the derivative...
    SimplePolynomial
    potentialDerivative( potentialApproximationOverAuxiliary );
    std::vector< double > const&
    derivativeVector( potentialDerivative.CoefficientVector() );
    for( unsigned int whichPower( 1 );
         whichPower < potentialApproximationPower;
         ++whichPower )
    {
      derivativeVector[ whichPower ] *= (double)( whichPower + 1 );
    }




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

} /* namespace VevaciousPlusPlus */
