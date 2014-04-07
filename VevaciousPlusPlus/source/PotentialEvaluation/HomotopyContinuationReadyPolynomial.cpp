/*
 * HomotopyContinuationReadyPolynomial.cpp
 *
 *  Created on: Apr 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  HomotopyContinuationReadyPolynomial::HomotopyContinuationReadyPolynomial() :
    HomotopyContinuationReadyPotential(),
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

} /* namespace VevaciousPlusPlus */
