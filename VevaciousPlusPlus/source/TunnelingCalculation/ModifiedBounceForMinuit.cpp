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
                                double const falseVacuumEvaporationTemperature,
                                unsigned int const undershootOvershootAttempts,
                             unsigned int const maximumMultipleOfLongestLength,
                                  double const initialFractionOfShortestLength,
                         unsigned int const energyConservingUndershootAttempts,
                                              double const minimumScaleSquared,
                                  double const shootingCloseEnoughThreshold ) :
    ROOT::Minuit2::FCNBase(),
    potentialFunction( potentialFunction ),
    numberOfFields( potentialFunction.NumberOfFieldVariables() ),
    referenceFieldIndex( 0 ),
    numberOfParameterizationFields( numberOfFields - 1 ),
    potentialApproximationPower( potentialApproximationPower ),
    falseVacuum( falseVacuum ),
    trueVacuum( trueVacuum ),
    fieldOriginPotential( potentialFunction(
                                         potentialFunction.FieldValuesOrigin(),
                                             0.0 ) ),
    falseVacuumPotential( potentialFunction( falseVacuum.FieldConfiguration(),
                                             0.0 ) ),
    trueVacuumPotential( potentialFunction( trueVacuum.FieldConfiguration(),
                                            0.0 ) ),
    falseVacuumEvaporationTemperature( falseVacuumEvaporationTemperature ),
    minimumScaleSquared( minimumScaleSquared ),
    tunnelingScaleSquared( std::max( minimumScaleSquared,
                            potentialFunction.ScaleSquaredRelevantToTunneling(
                                                                   falseVacuum,
                                                              trueVacuum ) ) ),
    shortestLength( 1.0 ),
    longestLength( 1.0 ),
    undershootOvershootAttempts( undershootOvershootAttempts ),
    energyConservingUndershootAttempts( energyConservingUndershootAttempts ),
    maximumMultipleOfLongestLength( maximumMultipleOfLongestLength ),
    initialFractionOfShortestLength( initialFractionOfShortestLength ),
    shootingThresholdSquared( shootingCloseEnoughThreshold
                              * shootingCloseEnoughThreshold )
  {
    // We pick the field that has the largest difference between the vacua as
    // the reference field that goes linearly with the auxiliary variable.
    std::vector< double > const&
    falseFieldConfiguration( falseVacuum.FieldConfiguration() );
    std::vector< double > const&
    trueFieldConfiguration( trueVacuum.FieldConfiguration() );
    double greatestFieldDifference( falseFieldConfiguration.front()
                                    - trueFieldConfiguration.front() );
    double currentFieldDifference( greatestFieldDifference );
    if( greatestFieldDifference < 0.0 )
    {
      greatestFieldDifference = -greatestFieldDifference;
    }
    for( size_t fieldIndex( 1 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentFieldDifference = ( falseFieldConfiguration[ fieldIndex ]
                                 - trueFieldConfiguration[ fieldIndex ] );
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


  // The bounce action is calculated by the following process:
  // 1) The path in field space from the false vacuum to the true vacuum is
  //    decoded from pathParameterization to give the field configuration f
  //    as a function of a path auxiliary variable p, giving f(p) and df/dp.
  //    See the comment above DecodePathParameters for the format of
  //    pathParameterization.
  // 2) The potential is fitted as a polynomial in p of degree
  //    potentialApproximationPower, giving V(p).
  // 3) If the potential difference between the vacua is small enough, the
  //    thin-wall approximation is tried, but only returned if the resulting
  //    bubble radius turns out to be large enough compared to the wall
  //    thickness.
  // 4) If the thin-wall approximation is not returned, the one-dimensional
  //    bubble equation of motion along p is integrated to get p as a
  //    function of the radial variable r, which is the spatial radius for
  //    three-dimensional actions, or the length of the four-dimensional
  //    Euclidean vector. This gets the correct bubble profile only if the
  //    transverse equations of motion in field space are also satisfied, or
  //    for a modified potential which has additional terms that raise the
  //    energy barrier on the side of the path that has a lower energy
  //    barrier. Hence this bubble profile gives an upper bound on the bounce
  //    action for the unmodified potential, as the modified potential (which
  //    has the same true vacuum and false vacuum and is exactly the same
  //    along the path) cannot have a lower bounce action.
  //    [The p equation of motion has a strange term giving a force
  //    proportional to (dp/dr)^2 which is not a friction term as the sign
  //    is not proportional to dp/dr, rather to d^2f/dp^2. This is because p
  //    is not linear in the distance in field space, so actually this term
  //    behaves a bit like extra kinetic energy because it rolled further
  //    "down the hill" from the true vacuum initially because there is more
  //    field space distance covered by a small change in p if the derivative
  //    is that way, or like extra friction because the approach to the
  //    false vacuum is actually longer in field space distance than in p.
  //    The alternative is to divide up the path into pathResolution segments
  //    and take the path as a series of straight lines at a constant
  //    "velocity" in field space, meaning that the weird pseudo-friction
  //    force in the bubble wall equation of motion in p disappears. This
  //    would require some contortions with linked lists (probably) of
  //    segments to get f(p(r)) and df/dp.]
  // 5) The bounce action along p(r) is numerically integrated and returned.
  double ModifiedBounceForMinuit::operator()(
                      std::vector< double > const& pathParameterization ) const
  {
    // Current outline (code below doesn't do most of this yet):
    // 1) Take pathParameterization as proposed dependence of field configuration
    //    f on auxiliary variable p. An extra power of p is given with a
    //    coefficient to bring f to the true vacuum at p = 1, as
    //    pathParameterization alone does not have a fixed endpoint at p = 1
    //    for different pathParameterization input.
    // 2) Evaluate the potential at ( potentialApproximationPower - 2 ) values
    //    for p between p = 0 and p = 1 inclusive, giving V(p) as a polynomial
    //    in p (though V(p) is an approximation, f(p) is exact by
    //    construction). (The number of points to evaluate V(p) would be
    //    ( potentialApproximationPower + 1 ) for potentialApproximationPower
    //    coefficients of non-zero powers of p plus a constant coefficient, but
    //    2 of the coefficients are fixed by the requirement that dV/dp is 0 at
    //    p = 0 and p = 1, and the constant coefficient by taking V(0) to be 0,
    //    which also is convenient for evaluating the action. Hence V(p) is
    //    actually a polynomial approximation of the potential difference from
    //    the false vacuum along the path given by f(p).)
    // 3a) Return early with a thin-wall result if appropriate.
    // 3b) Find the 1-dimensional bubble profile by solving the radial bubble
    //     equations of motion along the path from p = 0 to p = p_crit, where
    //     p_crit is the value of p for which the bubble is critical, and is
    //     found by the undershoot/overshoot method (0.0 < p_crit < 1.0). This
    //     probably will involve using the Boost library odeint.
    // 4) Combine f(p) (exact) and V(p) (approximate) with p(r) to numerically
    //    integrate the bounce action.
    // Note that this is only really the bounce action (within the validity of
    // the approximations of truncated expansions) if the equations of motion
    // perpendicular to f(p) are also satisfied. However, they would be
    // satisfied for a modified potential which increases the energy barrier on
    // the side of the path until the "force" on the path balances, and such a
    // modified potential still has the same true and false vacua, and its
    // decay width is given by this calculated bounce action. The unmodified
    // potential with the possibility for a path with a lower energy barrier
    // must necessarily have at least as large a decay width as the modified
    // potential, hence the calculated bounce action here is an upper bound on
    // the true bounce action, and minimizing this should lead one to the true
    // minimal bounce action.

    double const givenTemperature( pathParameterization.back() );
    bool const nonZeroTemperature( givenTemperature > 0.0 );

    std::vector< SimplePolynomial > fieldsAsPolynomials;
    double falseVacuumRelativePotential( falseVacuumPotential );
    double trueVacuumRelativePotential( trueVacuumPotential );
    DecodePathParameters( pathParameterization,
                          fieldsAsPolynomials,
                          falseVacuumRelativePotential,
                          trueVacuumRelativePotential );

    SimplePolynomial
    potentialApproximation( PotentialAlongPath( fieldsAsPolynomials,
                                                falseVacuumRelativePotential,
                                               trueVacuumRelativePotential ) );
    std::vector< SimplePolynomial > fieldDerivatives( numberOfFields );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      fieldDerivatives[ fieldIndex ]
      = fieldsAsPolynomials[ fieldIndex ].FirstDerivative();
    }

    // We return a thin-wall approximation if appropriate:
    double bounceAction( NAN );
    if( ThinWallAppropriate( ( falseVacuumRelativePotential
                               - trueVacuumRelativePotential ),
                             givenTemperature,
                             fieldDerivatives,
                             potentialApproximation,
                             bounceAction ) )
    {
      return bounceAction;
    }
    // If we didn't return bounceAction already, it means that the we go on to
    // try undershooting/overshooting.

    unsigned int dampingFactor( 3 );
    if( nonZeroTemperature )
    {
      dampingFactor = 2;
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
    double undershootAuxiliary( 0.0 );
    double overshootAuxiliary( 1.0 );
    // First we use conservation of energy to get a lower bound for the range
    // for the auxiliary variable which definitely undershoots.
    for( size_t undershootGuessStep( 0 );
         undershootGuessStep < energyConservingUndershootAttempts;
         ++undershootGuessStep )
    {
      if( potentialApproximation( 0.5 * ( undershootAuxiliary
                                          + overshootAuxiliary ) ) < 0.0 )
      {
        overshootAuxiliary
        = ( 0.5 * ( undershootAuxiliary + overshootAuxiliary ) );
      }
      else
      {
        undershootAuxiliary
        = ( 0.5 * ( undershootAuxiliary + overshootAuxiliary ) );
      }
    }
    // At this point we reset overshootAuxiliary for the integration including
    // damping.
    overshootAuxiliary = 1.0;
    double
    currentAuxiliary( 0.5 * ( undershootAuxiliary + overshootAuxiliary ) );
    size_t shootAttempts( 0 );
    double integrationRadius( longestLength );
    double
    initialIntegrationStep( initialFractionOfShortestLength * shortestLength );
    // size_t integrationSteps( 0 );
    std::vector< double > initialConditions( 2,
                                             0.0 );
    // The unit for the radial variable is 1/GeV.

    // This loop is broken out of if the shoot attempt seems to have been close
    // enough that the integration would take too long to find an overshoot or
    // undershoot, or that the shot was dead on.
    while( (++shootAttempts) <= undershootOvershootAttempts )
    {
      initialConditions[ 0 ] = currentAuxiliary;
      initialConditions[ 1 ] = 0.0;
      integrationRadius = longestLength;
      initialIntegrationStep
      = ( initialFractionOfShortestLength * shortestLength );
      bubbleProfiler.DoFirstStep( initialConditions,
                                  currentAuxiliary,
                                  initialIntegrationStep );
      odeintBubbleObserver.ResetValues( currentAuxiliary );
      // integrationSteps =
      boost::numeric::odeint::integrate( bubbleProfiler,
                                         initialConditions,
                                         initialIntegrationStep,
                                         integrationRadius,
                                         initialIntegrationStep,
                                         odeintBubbleObserver );
      odeintBubbleObserver.SortByRadialValue();

      while( WorthIntegratingFurther( odeintBubbleObserver,
                                      fieldsAsPolynomials ) )
      {
        integrationRadius *= 2.0;
        std::vector< BubbleRadialValueDescription >::reverse_iterator
        bubbleDescriptionIterator(
                           odeintBubbleObserver.BubbleDescription().rbegin() );
        initialConditions[ 0 ] = bubbleDescriptionIterator->auxiliaryValue;
        initialConditions[ 1 ] = bubbleDescriptionIterator->auxiliarySlope;
        initialIntegrationStep = bubbleDescriptionIterator->radialValue;
        initialIntegrationStep -= (++bubbleDescriptionIterator)->radialValue;
        // integrationSteps =
        boost::numeric::odeint::integrate( bubbleProfiler,
                                           initialConditions,
                                           ( 0.5 * integrationRadius ),
                                           integrationRadius,
                                           initialIntegrationStep,
                                           odeintBubbleObserver );
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

    double fieldGradient( 0.0 );
    double gradientDotGradient( 0.0 );

    bounceAction = 0.0;
    double previousIntegrand( 0.0 );
    // The bounce action density at r = 0 is 0 by merit of
    // r_0^dampingFactor = 0,
    // dp/dr = 0 at r = 0,
    // and potentialApproximation(p(0)) = 0 by construction.
    double currentIntegrand( 0.0 );
    double previousRadius( 0.0 );
    double currentRadius( 0.0 );
    double currentAuxiliary( 0.0 );
    double halfAuxiliarySlopeSquared( 0.0 );
    double fieldDerivative( 0.0 );
    double fieldDerivativeSquared( 0.0 );
    std::vector< BubbleRadialValueDescription > const&
    bubbleProfile( odeintBubbleObserver.BubbleDescription() );
    for( unsigned int radiusIndex( 1 );
         radiusIndex < bubbleProfile.size();
         ++radiusIndex )
    {
      previousRadius = currentRadius;
      previousIntegrand = currentIntegrand;
      currentRadius = bubbleProfile[ radiusIndex ].radialValue;
      currentAuxiliary = bubbleProfile[ radiusIndex ].auxiliaryValue;
      halfAuxiliarySlopeSquared = bubbleProfile[ radiusIndex ].auxiliarySlope;
      halfAuxiliarySlopeSquared *= ( 0.5 * halfAuxiliarySlopeSquared );
      fieldDerivativeSquared = 0.0;
      for( std::vector< SimplePolynomial >::const_iterator
           fieldDerivative( fieldDerivatives.begin() );
           fieldDerivative < fieldDerivatives.end();
           ++fieldDerivative )
      {
        fieldDerivative = (*fieldDerivative)( currentAuxiliary );
        fieldDerivativeSquared += ( fieldDerivative * fieldDerivative );
      }
      currentIntegrand
      = ( pow( currentRadius,
               dampingFactor )
          * ( halfAuxiliarySlopeSquared
              + potentialApproximation( currentAuxiliary ) ) );
      bounceAction += ( ( currentIntegrand + previousIntegrand )
                        * ( currentRadius - previousRadius ) );
      // A common factor of 1/2 is being left until after the loop.
    }
    // The common factor of 1/2 is combined with the solid angle of
    // 2 pi^2 (quantum) or 4 pi (thermal):
    if( nonZeroTemperature )
    {
      bounceAction *= ( 2.0 * M_PI );
    }
    else
    {
      bounceAction *= ( M_PI * M_PI );
    }

    if( nonZeroTemperature )
    {
      return ( ( bounceAction / givenTemperature ) + log( bounceAction ) );
    }
    else
    {
      return bounceAction;
    }
  }

  // This turns a flattened matrix of coefficients from pathParameterization
  // and fills fieldsAsPolynomials appropriately. The coefficients are taken
  // to be in the order
  // [ c_{1,0}, c_{1,1}, ..., c_{1, (referenceFieldIndex-1)},
  //           c_{1, (referenceFieldIndex+1)}, ..., c_{1,(numberOfFields-1)},
  //   c_{2,0}, c_{2,1}, ..., c_{2, (referenceFieldIndex-1)},
  //           c_{2, (referenceFieldIndex+1)}, ..., c_{1,(numberOfFields-1)},
  //   ...
  //   c_{p,0}, c_{d,1}, ..., c_{d, (referenceFieldIndex-1)},
  //           c_{d, (referenceFieldIndex+1)}, ..., c_{d,(numberOfFields-1)},
  //   temperature ],
  // where given field [j] is then the sum of c_{i,j} * p^i and d is the
  // greatest power given implicitly by pathParameterization. A final
  // coefficient for each polynomial is given so that the field takes its
  // value at the true vacuum for p = 1. Note that the given fields do not
  // map completely to the fields of potentialFunction: the field with index
  // referenceFieldIndex is skipped over by pathParameterization, as it is
  // set to be linear in a going from the false vacuum to the true vacuum. It
  // also puts the value of the potential (minus the value at the field
  // origin at zero temperature) into thermalFalseVacuumPotential and
  // thermalTrueVacuumPotential for the false and true vacua at the
  // temperature given by pathParameterization.back() if and only if it is
  // non-zero. Unfortunately it does not use "return value optimization" as
  // we cannot be sure that the user will compile with C++11 features
  // enabled.
  void ModifiedBounceForMinuit::DecodePathParameters(
                             std::vector< double > const& pathParameterization,
                          std::vector< SimplePolynomial >& fieldsAsPolynomials,
                                           double& thermalFalseVacuumPotential,
                                     double& thermalTrueVacuumPotential ) const
  {
    double const givenTemperature( pathParameterization.back() );
    bool const nonZeroTemperature( givenTemperature > 0.0 );

    std::vector< double >
    falseVacuumConfiguration( falseVacuum.FieldConfiguration() );
    std::vector< double >
    trueVacuumConfiguration( trueVacuum.FieldConfiguration() );

    if( nonZeroTemperature )
    {
      PotentialForMinuit potentialForMinuit( potentialFunction );
      MinuitManager thermalDsbFinder( potentialForMinuit );
      MinuitMinimum
      thermalTrueVacuum( thermalDsbFinder( trueVacuumConfiguration ) );
      trueVacuumConfiguration = thermalTrueVacuum.VariableValues();
      thermalTrueVacuumPotential = ( thermalTrueVacuum.FunctionValue()
                                     + potentialForMinuit.FunctionAtOrigin() );
      // We undo the offset from potentialFromMinuit.
      if( givenTemperature < falseVacuumEvaporationTemperature )
      {
        MinuitMinimum thermalFalseVacuum
        = thermalDsbFinder( falseVacuum.FieldConfiguration() );
        falseVacuumConfiguration = thermalFalseVacuum.VariableValues();
        thermalFalseVacuumPotential = ( thermalFalseVacuum.FunctionValue()
                                     + potentialForMinuit.FunctionAtOrigin() );
        // We undo the offset from potentialFromMinuit.
      }
      else
      {
        falseVacuumConfiguration = std::vector< double >( numberOfFields,
                                                          0.0 );
        // We use potentialForMinuit for a consistent offset to the value of
        // potentialFunction at the field origin at zero temperature.
        thermalFalseVacuumPotential
        = potentialFunction( falseVacuumConfiguration,
                             givenTemperature );
      }
    }

    size_t coefficientsPerField( ( ( pathParameterization.size() - 2 )
                                   / numberOfParameterizationFields ) + 3 );

    // At this point, we know that the path parameterization is a flattened
    // coefficientsPerField * numberOfParameterizationFields matrix.


    // Now we need to work out what we do!
    // Maybe this is best:
    // vector in field space from false to true is v
    // pathParameterization has sets of vectors in plane perpendicular to
    // reference field axis each p of those is projected onto plane
    // perpendicular to v by p' = p - ( p.v v / v^2 ), then has
    // v * index * numberOfParameterizationFields / pathParameterization.size()
    // added to it.
    // p'[i] = p + [ ( n * i / s ) - ( p.v * v^(-2) ) ] v


    // Once we have a set of nodes, we might want to project them from planes
    // perpendicular to v onto maybe hyperbolae? I dunno, something to maybe
    // allow the path to leave the endpoints in directions that initially
    // increase the distance to the other endpoint...

    fieldsAsPolynomials.resize( numberOfFields,
                                SimplePolynomial( coefficientsPerField ) );
    fieldsAsPolynomials[ referenceFieldIndex ].CoefficientVector().resize( 2,
                             falseVacuumConfiguration[ referenceFieldIndex ] );
    fieldsAsPolynomials[ referenceFieldIndex ].CoefficientVector().back()
    = ( trueVacuumConfiguration[ referenceFieldIndex ]
        - falseVacuumConfiguration[ referenceFieldIndex ] );
    // The last element of pathParameterization is a temperature, so should not
    // be used. The ( ( ( size - 2 ) / numberOfParameterizationFields ) + 3 ) integer is
    // to ensure that each SimplePolynomial has enough elements. The first
    // coefficient (the constant) is taken from falseVacuumConfiguration. The
    // last coefficient of each field is implicit, and is fixed by the
    // requirement that the field path goes to trueVacuumConfiguration for a=1.
    size_t fieldIndexFromSplines( 0 );
    for( size_t splineIndex( 0 );
         splineIndex < ( pathParameterization.size() - 1 );
         ++splineIndex )
    {
      fieldIndexFromSplines = ( splineIndex % numberOfParameterizationFields );
      if( fieldIndexFromSplines >= referenceFieldIndex )
      {
        ++fieldIndexFromSplines;
      }
      fieldsAsPolynomials[ fieldIndexFromSplines ].CoefficientVector()[
                         ( splineIndex / numberOfParameterizationFields ) + 1 ]
      = pathParameterization[ splineIndex ];
    }
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      std::vector< double >& coefficientVector(
                       fieldsAsPolynomials[ fieldIndex ].CoefficientVector() );
      coefficientVector.front() = falseVacuumConfiguration[ fieldIndex ];
      double coefficientSum( 0.0 );
      for( size_t coefficientIndex( 0 );
           coefficientIndex < ( coefficientVector.size() - 1 );
           ++coefficientIndex )
      {
        coefficientSum += coefficientVector[ coefficientIndex ];
      }
      coefficientVector.back() = ( trueVacuumConfiguration[ fieldIndex ]
                                   - coefficientSum );
      // This ensures that the field configuration for auxiliary variable = 1.0
      // goes to that of the true vacuum.
    }
  }

  // This returns a polynomial approximation of the potential along the path
  // given by splineCoefficients.
  SimplePolynomial ModifiedBounceForMinuit::PotentialAlongPath(
                    std::vector< SimplePolynomial > const& fieldsAsPolynomials,
                                             double const falseVacuumPotential,
                                       double const trueVacuumPotential ) const
  {
    std::vector< double > fieldConfiguration( numberOfFields );
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
      SetFieldConfiguration( fieldConfiguration,
                             fieldsAsPolynomials,
                             auxiliaryValue );
      potentialValues( whichStep ) = ( potentialFunction( fieldConfiguration )
                                       - falseVacuumPotential );

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
    = ( trueVacuumPotential - falseVacuumPotential );
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
                            2,
                            1 );
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
    potentialVector.back() = ( finalCoefficientTimesMaxPower
                               / (double)potentialApproximationPower );
    return potentialApproximation;
  }

  // This puts the bounce action in the thin-wall approximation into
  // bounceAction and returns true if the thin-wall approximation is
  // appropriate. Otherwise, bounceAction is left alone and false is
  // returned.
  bool
  ModifiedBounceForMinuit::ThinWallAppropriate( double potentialDifference,
                                                double const givenTemperature,
                       std::vector< SimplePolynomial > const& fieldDerivatives,
                                SimplePolynomial const& potentialApproximation,
                                                double& bounceAction ) const
  {
    if( potentialDifference < 0.0 )
    {
      potentialDifference = -potentialDifference;
    }
    if( potentialDifference
        > ( 1.0E-3 * tunnelingScaleSquared * tunnelingScaleSquared ) )
    {
      return false;
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
      potentialValue = potentialApproximation( thinWallAuxiliary );
      derivativeSquared = 0.0;
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        fieldDerivative = fieldDerivatives[ fieldIndex ]( thinWallAuxiliary );
        derivativeSquared += ( fieldDerivative * fieldDerivative );
      }
      integratedAction += sqrt( derivativeSquared * potentialValue );
      wallThickness += ( thinWallStep / sqrt( potentialValue ) );
    }
    integratedAction *= ( M_SQRT2 * thinWallStep );
    wallThickness
    *= fieldDerivatives[ referenceFieldIndex ].CoefficientVector().front();
    // This last *= accounts for the change in integration variable to the
    // field linear in a.


    // We return the thin-wall approximation of the bounce action only if the
    // bubble radius is sufficiently large compared to the tunneling scale,
    // which is a good indicator for the validity of the approximation.
    if( ( integratedAction * oneOverPotentialDifference )
        < ( 1.0E+2 * wallThickness ) )
    {
      return false;
    }
    if( givenTemperature > 0.0 )
    {
      bounceAction = ( ( 2.0 * M_PI * integratedAction * integratedAction
                       * integratedAction * oneOverPotentialDifference
                       * oneOverPotentialDifference ) / givenTemperature );
      // This is at least 10^4 * 2 pi * ( S_1 / ( Q^2 T ) ), so S_3(T)/T is
      // likely to be at least 10^3, leading to a totally negligible
      // tunneling probability...
    }
    else
    {
      bounceAction = ( -13.5 * M_PI * M_PI * integratedAction
                       * integratedAction * integratedAction * integratedAction
                       * oneOverPotentialDifference
                       * oneOverPotentialDifference
                       * oneOverPotentialDifference );
      // This is at least 3^4 * 10^4 * 13.5 pi^2 * ( S_1 / ( Q^3 ) ), so
      // S_4 is likely to be at least 10^6, leading to a totally negligible
      // tunneling probability...
    }
    return true;
  }

  // This returns true if neither overshooting nor undershooting definitely
  // happened and the distance in field space to the top of the false vacuum
  // hill is still larger than shootingCloseEnoughThreshold times the
  // distance between the initial field configuration and the false vacuum
  // field configuration.
  bool ModifiedBounceForMinuit::WorthIntegratingFurther(
                              OdeintBubbleObserver const& odeintBubbleObserver,
              std::vector< SimplePolynomial >const& fieldsAsPolynomials ) const
  {
    if( odeintBubbleObserver.DefinitelyUndershot()
        ||
        odeintBubbleObserver.DefinitelyOvershot() )
    {
      return false;
    }
    double totalDistanceSquared( 0.0 );
    double nearDistanceSquared( 0.0 );
    double fieldDifference( 0.0 );
    double falseVacuumField( 0.0 );
    double const nearAuxiliary(
              odeintBubbleObserver.BubbleDescription().back().auxiliaryValue );
    double const startAuxiliary(
             odeintBubbleObserver.BubbleDescription().front().auxiliaryValue );
    for( std::vector< SimplePolynomial >::const_iterator
         fieldValue( fieldsAsPolynomials.begin() );
         fieldValue < fieldsAsPolynomials.end();
         ++fieldValue )
    {
      falseVacuumField = fieldValue->CoefficientVector().front();
      fieldDifference = ( (*fieldValue)( nearAuxiliary ) - falseVacuumField );
      nearDistanceSquared += ( fieldDifference * fieldDifference );
      fieldDifference = ( (*fieldValue)( startAuxiliary ) - falseVacuumField );
      totalDistanceSquared += ( fieldDifference * fieldDifference );
    }
    return ( nearDistanceSquared
             > ( shootingThresholdSquared * totalDistanceSquared ) );
  }

} /* namespace VevaciousPlusPlus */
