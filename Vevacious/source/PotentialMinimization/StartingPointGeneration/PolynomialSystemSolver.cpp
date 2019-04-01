/*
 * PolynomialSystemSolver.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/StartingPointGeneration/PolynomialSystemSolver.hpp"

namespace VevaciousPlusPlus
{

  PolynomialSystemSolver::PolynomialSystemSolver()
  {
    // This does nothing.
  }

  PolynomialSystemSolver::~PolynomialSystemSolver()
  {
    // This does nothing.
  }


  // This goes through each of the constraints in systemToSolve, stepping
  // resolutionSize either side of givenSolution in the field appropriate to
  // the constraint, and returns false if any of the constraints do not change
  // sign in stepping from one side of givenSolution to the other in any of the
  // fields.
  bool PolynomialSystemSolver::IsValidSolution(
                                           std::vector< double > givenSolution,
                      std::vector< PolynomialConstraint > const& systemToSolve,
                                                double const resolutionSize )
  {
    size_t const numberOfElements( givenSolution.size() );
    double originalValue( 0.0 );
    double negativePartialSlope( 0.0 );
    double positivePartialSlope( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfElements;
         ++fieldIndex )
    {
      originalValue = givenSolution[ fieldIndex ];
      givenSolution[ fieldIndex ] = ( originalValue - resolutionSize );
      negativePartialSlope = PartialSlope( systemToSolve[ fieldIndex ],
                                           givenSolution );
      givenSolution[ fieldIndex ] = ( originalValue + resolutionSize );
      positivePartialSlope = PartialSlope( systemToSolve[ fieldIndex ],
                                           givenSolution );

      if( ( ( positivePartialSlope > 0.0 )
            &&
            ( negativePartialSlope > 0.0 ) )
          ||
          ( ( positivePartialSlope < 0.0 )
            &&
            ( negativePartialSlope < 0.0 ) ) )
      {
        return false;
      }
      givenSolution[ fieldIndex ] = originalValue;
    }
    return true;
  }

  // This returns the value of fieldConstraint for the field values given in
  // fieldConfiguration.
  double PolynomialSystemSolver::PartialSlope(
                                   PolynomialConstraint const& fieldConstraint,
                                     std::vector< double > fieldConfiguration )
  {
    double partialSlope( 0.0 );
    double polynomialValue( 0.0 );
    for( std::vector< FactorWithPowers >::const_iterator
         factorWithPowers( fieldConstraint.begin() );
         factorWithPowers != fieldConstraint.end();
         ++factorWithPowers )
    {
      polynomialValue = factorWithPowers->first;
      for( size_t fieldIndex( 0 );
           fieldIndex < factorWithPowers->second.size();
           ++fieldIndex )
      {
        if( factorWithPowers->second[ fieldIndex ] > 0 )
        {
          switch( factorWithPowers->second[ fieldIndex ] )
          {
            case 1:
              polynomialValue *= fieldConfiguration[ fieldIndex ];
              break;
            case 2:
              polynomialValue *= ( fieldConfiguration[ fieldIndex ]
                                   * fieldConfiguration[ fieldIndex ] );
              break;
            case 3:
              polynomialValue *= ( fieldConfiguration[ fieldIndex ]
                                   * fieldConfiguration[ fieldIndex ]
                                   * fieldConfiguration[ fieldIndex ] );
              break;
            default:
              for( size_t powerCount( 0 );
                   powerCount < factorWithPowers->second[ fieldIndex ];
                   ++powerCount )
              {
                polynomialValue *= fieldConfiguration[ fieldIndex ];
              }
              break;
          }
        }
      }
      partialSlope += polynomialValue;
    }
    return partialSlope;
  }

}
