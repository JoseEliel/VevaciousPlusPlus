/*
 * HomotopyContinuationReadyPolynomial.hpp
 *
 *  Created on: Apr 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONREADYPOLYNOMIAL_HPP_
#define HOMOTOPYCONTINUATIONREADYPOLYNOMIAL_HPP_

#include "../../StandardIncludes.hpp"
#include "HomotopyContinuationReadyPotential.hpp"
#include "../PolynomialTerm.hpp"
#include "../PolynomialSum.hpp"
#include "../ProductOfPolynomialSums.hpp"

namespace VevaciousPlusPlus
{

  class HomotopyContinuationReadyPolynomial :
                                      public HomotopyContinuationReadyPotential
  {
  public:
    HomotopyContinuationReadyPolynomial( SlhaManager& slhaManager );
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

    // This evaluates the target system and places the values in
    // destinationVector.
    virtual void HomotopyContinuationSystemValues(
                   std::vector< std::complex< double > > solutionConfiguration,
                    std::vector< std::complex< double > >& destinationVector );

    // This evaluates the derivatives of the target system and places the
    // values in destinationMatrix.
    virtual void HomotopyContinuationSystemGradients(
                   std::vector< std::complex< double > > solutionConfiguration,
                                 std::vector< std::vector< std::complex< double
                                                    > > >& destinationMatrix );


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

  // This evaluates the target system and places the values in
  // destinationVector.
  inline void
  HomotopyContinuationReadyPolynomial::HomotopyContinuationSystemValues(
                   std::vector< std::complex< double > > solutionConfiguration,
                     std::vector< std::complex< double > >& destinationVector )
  {
    destinationVector.resize( targetPolynomialGradient.size() );
    for( unsigned int whichIndex( 0 );
         whichIndex < targetPolynomialGradient.size();
         ++whichIndex )
    {
      destinationVector[ whichIndex ]
      = targetPolynomialGradient[ whichIndex ]( solutionConfiguration );
    }
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "HomotopyContinuationReadyPolynomial::HomotopyContinuationSystemValues("
    << " {";
    for( unsigned int whichIndex( 0 );
         whichIndex < solutionConfiguration.size();
         ++whichIndex )
    {
      std::cout << " ( " << solutionConfiguration[ whichIndex ].real()
       << ", " << solutionConfiguration[ whichIndex ].imag() << " )";
    }
    std::cout << " }, ... ) set destinationVector to be { ";
    for( unsigned int whichIndex( 0 );
         whichIndex < destinationVector.size();
         ++whichIndex )
    {
      std::cout << " ( " << destinationVector[ whichIndex ].real()
       << ", " << destinationVector[ whichIndex ].imag() << " )";
    }
    std::cout << " }.";
    std::cout << std::endl;
    std::cout << std::endl;*/
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

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "HomotopyContinuationReadyPolynomial::PreparePolynomialGradient()"
    <<  " called. targetPolynomialGradient.size() set to "
    << targetPolynomialGradient.size()
    <<  ", polynomialForGradient.size() = " << polynomialForGradient.size();
    std::cout << std::endl;*/

    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      for( std::vector< PolynomialTerm >::const_iterator
           whichTerm( polynomialForGradient.begin() );
           whichTerm < polynomialForGradient.end();
           ++whichTerm )
      {
        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "trying to differentiate [" << whichTerm->AsDebuggingString()
        << "] with respect to fieldNames[ " << fieldIndex << " ] = \""
        << fieldNames[ fieldIndex ] << "\".";
        std::cout << std::endl;*/

        if( whichTerm->NonZeroDerivative( fieldIndex ) )
        {
          targetPolynomialGradient[ fieldIndex ].PolynomialTerms().push_back(
                                  whichTerm->PartialDerivative( fieldIndex ) );
        }
      }
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "targetPolynomialGradient =" << std::endl;
    for( std::vector< PolynomialSum >::iterator
         whichTadpole( targetPolynomialGradient.begin() );
         whichTadpole < targetPolynomialGradient.end();
         ++whichTadpole )
    {
      std::cout << std::endl << whichTadpole->AsDebuggingString();
    }
    std::cout << std::endl;
    std::cout << std::endl;*/
  }

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONREADYPOLYNOMIAL_HPP_ */
