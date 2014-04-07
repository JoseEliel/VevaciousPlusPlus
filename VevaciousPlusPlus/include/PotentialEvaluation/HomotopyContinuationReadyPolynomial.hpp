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
#include "ProductOfPolynomialSums.hpp"

namespace VevaciousPlusPlus
{

  class HomotopyContinuationReadyPolynomial :
                                      public HomotopyContinuationReadyPotential
  {
  public:
    HomotopyContinuationReadyPolynomial();
    virtual
    ~HomotopyContinuationReadyPolynomial();


    // This fills targetPolynomialGradient, homotopyContinuationStartSystem,
    // targetPolynomialHessian, and startPolynomialHessian appropriately.
    virtual void
    PreparePolynomialHomotopyContinuation();

    std::vector< PolynomialSum > const& TargetPolynomialGradient() const
    { return targetPolynomialGradient; }

    std::vector< std::vector< PolynomialSum > > const&
    TargetPolynomialHessian() const{ return targetPolynomialHessian; }

    std::vector< ProductOfPolynomialSums > const&
    HomotopyContinuationStartSystem() const
    { return homotopyContinuationStartSystem; }

    std::vector< std::vector< std::complex< double > > > const&
    HomotopyContinuationStartValues() const
    { return homotopyContinuationStartValues; }

    std::vector< std::vector< ProductOfPolynomialSums > > const&
    HomotopyContinuationStartHessian() const
    { return startPolynomialHessian; }


  protected:
    PolynomialSum homotopyContinuationPotentialPolynomial;
    std::vector< PolynomialSum > targetPolynomialGradient;
    std::vector< ProductOfPolynomialSums > homotopyContinuationStartSystem;
    std::vector< std::vector< std::complex< double > > >
    homotopyContinuationStartValues;
    std::vector< std::vector< std::complex< double > > >
    homotopyContinuationValidSolutions;
    std::vector< std::vector< PolynomialSum > > targetPolynomialHessian;
    std::vector< std::vector< ProductOfPolynomialSums > >
    startPolynomialHessian;

    // This should prepare homotopyContinuationPotentialPolynomial
    // appropriately.
    virtual void PrepareHomotopyContinuationPotentialPolynomial() = 0;

    // This fills targetPolynomialGradient from
    // homotopyContinuationPotentialPolynomial.
    virtual void PreparePolynomialGradient();

    // This fills targetPolynomialHessian from targetPolynomialGradient.
    virtual void PreparePolynomialHessian();

    // This should prepare homotopyContinuationStartSystem and
    // startPolynomialHessian appropriately.
    virtual void PrepareHomotopyContinuationStartSystem() = 0;

    // This should prepare homotopyContinuationStartValues to be all the
    // solutions of homotopyContinuationStartSystem.
    virtual void PrepareHomotopyContinuationStartValues() = 0;
  };




  // This fills targetPolynomialGradient, homotopyContinuationStartSystem,
  // targetPolynomialHessian, and startPolynomialHessian appropriately.
  inline void
  HomotopyContinuationReadyPolynomial::PreparePolynomialHomotopyContinuation()
  {
    PrepareHomotopyContinuationPotentialPolynomial();
    PreparePolynomialGradient();
    PreparePolynomialHessian();
    PrepareHomotopyContinuationStartSystem();
    PrepareHomotopyContinuationStartValues();
  }

  // This fills targetPolynomialGradient from
  // homotopyContinuationPotentialPolynomial.
  inline void
  HomotopyContinuationReadyPolynomial::PreparePolynomialGradient()
  {
    targetPolynomialGradient.assign( numberOfFields,
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
          targetPolynomialGradient[ fieldIndex ].PolynomialTerms().push_back(
                                  whichTerm->PartialDerivative( fieldIndex ) );
        }
      }
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONREADYPOLYNOMIAL_HPP_ */
