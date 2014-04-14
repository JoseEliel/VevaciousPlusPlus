/*
 * HomotopyContinuationReadyPolynomial.cpp
 *
 *  Created on: Apr 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  HomotopyContinuationReadyPolynomial::HomotopyContinuationReadyPolynomial(
                                                   SlhaManager& slhaManager ) :
    HomotopyContinuationReadyPotential( slhaManager ),
    homotopyContinuationPotentialPolynomial(),
    targetPolynomialGradient(),
    homotopyContinuationStartSystem(),
    homotopyContinuationStartValues(),
    homotopyContinuationValidSolutions(),
    targetPolynomialHessian(),
    startPolynomialHessian()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "HomotopyContinuationReadyPolynomial::"
    << "HomotopyContinuationReadyPolynomial()";
    std::cout << std::endl;/**/
  }

  HomotopyContinuationReadyPolynomial::~HomotopyContinuationReadyPolynomial()
  {
    // This does nothing.
  }


  // This fills targetPolynomialHessian from targetPolynomialGradient.
  void
  HomotopyContinuationReadyPolynomial::PreparePolynomialHessian()
  {
    targetPolynomialHessian.resize( numberOfFields,
                              std::vector< PolynomialSum >( numberOfFields ) );
    for( unsigned int gradientIndex( 0 );
         gradientIndex < numberOfFields;
         ++gradientIndex )
    {
      std::vector< PolynomialTerm > const& gradientVector(
                 targetPolynomialGradient[ gradientIndex ].PolynomialTerms() );
      for( unsigned int fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        for( std::vector< PolynomialTerm >::const_iterator
             whichTerm( gradientVector.begin() );
             whichTerm < gradientVector.end();
             ++whichTerm )
        {
          if( whichTerm->NonZeroDerivative( fieldIndex ) )
          {
            targetPolynomialHessian[ gradientIndex ][ fieldIndex
                                                             ].PolynomialTerms(
                     ).push_back( whichTerm->PartialDerivative( fieldIndex ) );
          }
        }
      }
    }
  }

  // This evaluates the derivatives of the target system and places the
  // values in destinationMatrix.
  void
  HomotopyContinuationReadyPolynomial::HomotopyContinuationSystemGradients(
                   std::vector< std::complex< double > > solutionConfiguration,
                                 std::vector< std::vector< std::complex< double
                                                     > > >& destinationMatrix )
  {
    destinationMatrix.resize( targetPolynomialHessian.size() );
    for( unsigned int constraintIndex( 0 );
         constraintIndex < targetPolynomialHessian.size();
         ++constraintIndex )
    {
      destinationMatrix[ constraintIndex ].resize(
                           targetPolynomialHessian[ constraintIndex ].size() );
      for( unsigned int variableIndex( 0 );
           variableIndex < targetPolynomialHessian[ constraintIndex ].size();
           ++variableIndex )
      {
        destinationMatrix[ constraintIndex ][ variableIndex ]
        = targetPolynomialHessian[ constraintIndex ][ variableIndex ](
                                                       solutionConfiguration );
      }
    }
  }

} /* namespace VevaciousPlusPlus */
