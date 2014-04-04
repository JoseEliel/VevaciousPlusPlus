/*
 * HomotopyContinuationReadyPolynomial.hpp
 *
 *  Created on: Apr 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONREADYPOLYNOMIAL_HPP_
#define HOMOTOPYCONTINUATIONREADYPOLYNOMIAL_HPP_

#include "../StandardIncludes.hpp"
#include "HomotopyContinuationReadyPotential.hpp"
#include "PolynomialTerm.hpp"
#include "PolynomialSum.hpp"

namespace VevaciousPlusPlus
{

  class HomotopyContinuationReadyPolynomial :
                                      public HomotopyContinuationReadyPotential
  {
  public:
    HomotopyContinuationReadyPolynomial();
    virtual
    ~HomotopyContinuationReadyPolynomial();


    // This fills polynomialGradient, homotopyContinuationStartSystem, and
    // polynomialHessian appropriately.
    virtual void
    PreparePolynomialHomotopyContinuation();

    std::vector< PolynomialSum > const& PolynomialGradient() const
    { return polynomialGradient; }

    std::vector< PolynomialSum > const&
    HomotopyContinuationStartSystem() const
    { return homotopyContinuationStartSystem; }

    std::vector< std::vector< std::complex< double > > > const&
    HomotopyContinuationStartValues() const
    { return homotopyContinuationStartValues; }

    std::vector< std::vector< PolynomialSum > > const&
    PolynomialHessian() const{ return polynomialHessian; }


  protected:
    PolynomialSum homotopyContinuationPotentialPolynomial;
    std::vector< PolynomialSum > polynomialGradient;
    std::vector< ProductOfPolynomialSums > homotopyContinuationStartSystem;
    std::vector< std::vector< std::complex< double > > >
    homotopyContinuationStartValues;
    std::vector< std::vector< ProductOfPolynomialSums > > polynomialHessian;

    // This should prepare homotopyContinuationPotentialPolynomial
    // appropriately.
    virtual void PrepareHomotopyContinuationPotentialPolynomial() = 0;

    // This fills polynomialGradient from
    // homotopyContinuationPotentialPolynomial.
    virtual void PrepareHomotopyContinuationGradient();

    // This fills polynomialHessian from polynomialGradient.
    virtual void PrepareHomotopyContinuationHessian();

    // This should prepare homotopyContinuationStartSystem appropriately.
    virtual void PrepareHomotopyContinuationStartSystem() = 0;

    // This should prepare homotopyContinuationStartValues to be all the
    // solutions of homotopyContinuationStartSystem.
    virtual void PrepareHomotopyContinuationStartValues() = 0;
  };




  // This fills polynomialGradient, homotopyContinuationStartSystem, and
  // polynomialHessian appropriately.
  inline void
  HomotopyContinuationReadyPolynomial::PreparePolynomialHomotopyContinuation()
  {
    PrepareHomotopyContinuationPotentialPolynomial();
    PrepareHomotopyContinuationGradient();
    PrepareHomotopyContinuationHessian();
    PrepareHomotopyContinuationStartSystem();
    PrepareHomotopyContinuationStartValues();
  }

  // This fills polynomialGradient from
  // homotopyContinuationPotentialPolynomial.
  inline void
  HomotopyContinuationReadyPolynomial::PrepareHomotopyContinuationGradient()
  {
    polynomialGradient.assign( numberOfFields,
                               PolynomialSum() );
    std::vector< PolynomialTerm > const& polynomialForGradient(
                   homotopyContinuationPotentialPolynomial.PolynomialTerms() );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      for( std::vector< PolynomialTerm >::const_iterator
           whichTerm( polynomialForGradient.begin() );
           whichTerm < polynomialForGradient.end();
           ++whichTerm )
      {
        if( whichTerm->NonZeroDerivative( fieldIndex ) )
        {
          polynomialGradient[ fieldIndex ].PolynomialTerms().push_back(
                                  whichTerm->PartialDerivative( fieldIndex ) );
        }
      }
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONREADYPOLYNOMIAL_HPP_ */
