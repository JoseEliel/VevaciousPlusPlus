/*
 * FieldPolynomialsWithScale.cpp
 *
 *  Created on: Apr 16, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/HomotopyContinuation/FieldPolynomialsWithScale.hpp"

namespace VevaciousPlusPlus
{

  FieldPolynomialsWithScale::FieldPolynomialsWithScale(
                                      PolynomialSum const& potentialPolynomial,
                                                   size_t const numberOfFields,
                                      SlhaUpdatePropagator& previousPropagator,
                            std::vector< size_t > const& fieldsAssumedPositive,
                            std::vector< size_t > const& fieldsAssumedNegative,
                                bool const treeLevelMinimaOnlyAsValidSolutions,
                            double const assumedPositiveOrNegativeTolerance ) :
    PolynomialGradientTargetSystem( potentialPolynomial,
                                    ( numberOfFields + 1 ),
                                    previousPropagator,
                                    fieldsAssumedPositive,
                                    fieldsAssumedNegative,
                                    treeLevelMinimaOnlyAsValidSolutions,
                                    assumedPositiveOrNegativeTolerance )
  {
    this->numberOfFields = numberOfFields;
  }

  FieldPolynomialsWithScale::~FieldPolynomialsWithScale()
  {
    // This does nothing.
  }


  // This calls the base PolynomialGradientTargetSystem version then replaces
  // targetSystem[ numberOfFields ] (which should be the scale constraint)
  // with a polynomial approximation of the constraint l = ln( Q ) where l is
  // the variable with index numberOfFields and Q is the scale taken as the
  // square root of the sum of the squares of the field values. Hence it
  // inserts a polynomial approximation of
  // exp( 2 * l ) = minimumScale^2 + f[0]^2 + f[1]^2 + ...
  // for the field configuration f[]. The minimum renormalization scale should
  // be given by lowerEndOfStartValues, and upperEndOfStartValues should give
  // the maximum renormalization scale.
  void FieldPolynomialsWithScale::UpdateSelfForNewSlha(
                                               SlhaManager const& slhaManager )
  {
    minimumScale = slhaManager.LowestBlockScale();
    maximumScale = slhaManager.HighestBlockScale();
    PolynomialGradientTargetSystem::UpdateSelfForNewSlha( slhaManager );

    // We approximate exp( 2 * l ) as
    // ( L - l )^(-1) * ( L + b * l + ( c / L ) * l^2 )
    // matched at l = 0, L/3, 2L/3, where L = ln(maximumScale)
    // b = 1.5 - 4.0 * exp( 2.0 / 3.0 ) + 1.5 * exp( 4.0 / 3.0 )
    // c = 4.5 * ( 4.0 * exp( 2.0 / 3.0 ) - 1.0 * exp( 4.0 / 3.0 ) - 3.0 )
    // Hence the constraint is
    // L + b * l + ( c / L ) * l^2
    // + l * f[1]^2 - L * f[1]^2 + l * f[2]^2 - L * f[2]^2 + ...
    double const logarithmOfMaximumScale( log( maximumScale ) );
    PolynomialSum scaleConstraint;
    PolynomialTerm freshTerm;
    freshTerm.MultiplyBy( logarithmOfMaximumScale );
    scaleConstraint.PolynomialTerms().push_back( freshTerm );
    freshTerm = PolynomialTerm();
    freshTerm.MultiplyBy( 1.5 - 4.0 * exp( 2.0 / 3.0 )
                              + 1.5 * exp( 4.0 / 3.0 ) );
    freshTerm.RaiseFieldPower( numberOfFields,
                               1 );
    scaleConstraint.PolynomialTerms().push_back( freshTerm );
    freshTerm = PolynomialTerm();
    freshTerm.MultiplyBy( ( 4.5 * ( 4.0 * exp( 2.0 / 3.0 )
                                    - exp( 4.0 / 3.0 ) - 3.0 ) )
                          / logarithmOfMaximumScale );
    freshTerm.RaiseFieldPower( numberOfFields,
                               2 );
    scaleConstraint.PolynomialTerms().push_back( freshTerm );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      freshTerm = PolynomialTerm();
      freshTerm.RaiseFieldPower( fieldIndex,
                                 2 );
      freshTerm.RaiseFieldPower( numberOfFields,
                                 1 );
      scaleConstraint.PolynomialTerms().push_back( freshTerm );
      freshTerm = PolynomialTerm();
      freshTerm.RaiseFieldPower( fieldIndex,
                                 2 );
      freshTerm.MultiplyBy( -logarithmOfMaximumScale );
      scaleConstraint.PolynomialTerms().push_back( freshTerm );
    }
    SetStartSystem( numberOfFields,
                    log( minimumScale ),
                    logarithmOfMaximumScale );
  }

} /* namespace VevaciousPlusPlus */
