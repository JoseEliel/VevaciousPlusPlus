/*
 * PolynomialGradientTargetSystem.hpp
 *
 *  Created on: Apr 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALGRADIENTTARGETSYSTEM_HPP_
#define POLYNOMIALGRADIENTTARGETSYSTEM_HPP_

#include "../../StandardIncludes.hpp"
#include "HomotopyContinuationTargetSystem.hpp"
#include "../../PotentialEvaluation/PolynomialTerm.hpp"
#include "../../PotentialEvaluation/PolynomialSum.hpp"
#include "../../PotentialEvaluation/ProductOfPolynomialSums.hpp"
#include "../../PotentialEvaluation/SumOfProductOfPolynomialSums.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialGradientTargetSystem :
                                        public HomotopyContinuationTargetSystem
  {
  public:
    PolynomialGradientTargetSystem( unsigned int const numberOfVariables );
    virtual
    ~PolynomialGradientTargetSystem();


    // This fills targetSystem, startSystem, targetHessian, and startHessian
    // appropriately. The starting system is set with
    // SetStartSystem( variableIndex,
    //                      lowerEndOfStartValues,
    //                      upperEndOfStartValues )
    // for variableIndex from 0 to ( numberOfVariables - 1 ). SetStartSystem
    // is public, so special ranges for particular variables can be given to
    // over-write the values given by this function.
    virtual void
    PrepareForHomotopyContinuation( PolynomialSum const& potentialPolynomial,
                                    double const lowerEndOfStartValues,
                                    double const upperEndOfStartValues );

    std::vector< PolynomialSum > const& TargetSystem() const
    { return targetSystem; }

    std::vector< PolynomialSum >& TargetSystem(){ return targetSystem; }

    std::vector< std::vector< PolynomialSum > > const& TargetHessian() const
    { return targetHessian; }

    std::vector< ProductOfPolynomialSums > const& StartSystem() const
    { return startSystem; }

    std::vector< std::vector< std::complex< double > > > const&
    StartValues() const{ return startValues; }

    std::vector< std::vector< ProductOfPolynomialSums > > const&
    StartHessian() const{ return startHessian; }

    // This evaluates the target system and places the values in
    // destinationVector.
    virtual void HomotopyContinuationSystemValues(
                   std::vector< std::complex< double > > solutionConfiguration,
                    std::vector< std::complex< double > >& destinationVector );

    // This evaluates the derivatives of the target system and places the
    // values in destinationMatrix.
    virtual void HomotopyContinuationSystemGradients(
                   std::vector< std::complex< double > > solutionConfiguration,
     std::vector< std::vector< std::complex< double > > >& destinationMatrix );

    // This sets startSystem[ variableIndex ] to be a polynomial of degree
    // highestDegreeOfTargetSystem with solutions given by
    // startValues[ variableIndex ] which is also set.
    // After that, startHessian[ variableIndex ] is also set accordingly.
    // It sets startSystem[ variableIndex ] to be of the form
    // x * ( x^2 - y_0^2 ) * ( x^2 - y_1^2 ) * ...
    // where x is the variable with index variableIndex and y_n is
    // startValues[ variableIndex ][ n ].
    // If highestDegreeOfTargetSystem is even,
    // startValues[ variableIndex ][ n-1 ] is set to n * q_n
    // for n going from 1 to d, where
    // d = highestDegreeOfTargetSystem / 2, and q_n is
    // lowerEndOfStartValues if d = 1 or
    // (lowerEndOfStartValues^([d-n]/[d-1])
    //  * upperEndOfStartValues^([n-1]/[d-1])) if d > 1.
    // If highestDegreeOfTargetSystem is odd, startValues[ variableIndex ][ n ]
    // is set by the above formula and startValues[ variableIndex ][ 0 ] is set
    // to 0.
    // So, if the highest degree polynomial of targetSystem is 3, d = 1, so the
    // starting values are set to
    // { 0, lowerEndOfStartValues }.
    // If highestDegreeOfTargetSystem is 4, d = 2, so:
    // { lowerEndOfStartValues, upperEndOfStartValues }.
    // If highestDegreeOfTargetSystem is 5, d = 2, so:
    // { 0, lowerEndOfStartValues, upperEndOfStartValues }.
    // If highestDegreeOfTargetSystem is 7, d = 3, so:
    // { 0, lowerEndOfStartValues,
    // ( lowerEndOfStartValues^(1/2) * upperEndOfStartValues^(1/2) ),
    // upperEndOfStartValues }.
    // If highestDegreeOfTargetSystem is 9, d = 4, so:
    // { 0, lowerEndOfStartValues,
    // ( lowerEndOfStartValues^(2/3) * upperEndOfStartValues^(1/3) ),
    // ( lowerEndOfStartValues^(1/3) * upperEndOfStartValues^(2/3) ),
    // upperEndOfStartValues }.
    void SetStartSystem( unsigned int const variableIndex,
                         double const lowerEndOfStartValues,
                         double const upperEndOfStartValues );


  protected:
    unsigned int const numberOfVariables;
    std::vector< PolynomialSum > targetSystem;
    unsigned int highestDegreeOfTargetSystem;
    std::vector< ProductOfPolynomialSums > startSystem;
    std::vector< std::vector< std::complex< double > > > startValues;
    std::vector< std::vector< std::complex< double > > > validSolutions;
    std::vector< std::vector< PolynomialSum > > targetHessian;
    std::vector< std::vector< SumOfProductOfPolynomialSums > > startHessian;

    // This fills targetSystem from potentialPolynomial.
    virtual void
    PreparePolynomialGradient( PolynomialSum const& potentialPolynomial );

    // This fills targetHessian from targetSystem.
    virtual void PreparePolynomialHessian();

    // This fills startHessian[ variableIndex ] from
    // startValues[ variableIndex ].
    void SetStartHessian( unsigned int const variableIndex );
  };




  // This fills targetSystem, startSystem, targetHessian, and startHessian
  // appropriately. The starting system is set with
  // SetStartSystem( variableIndex,
  //                 lowerEndOfStartValues,
  //                 upperEndOfStartValues )
  // for variableIndex from 0 to ( numberOfVariables - 1 ).
  inline void PolynomialGradientTargetSystem::PrepareForHomotopyContinuation(
                                      PolynomialSum const& potentialPolynomial,
                                            double const lowerEndOfStartValues,
                                          double const upperEndOfStartValues )
  {
    highestDegreeOfTargetSystem
    = ( potentialPolynomial.HighestFieldPower() - 1 );
    PreparePolynomialGradient( potentialPolynomial );
    PreparePolynomialHessian();
    startSystem.resize( numberOfVariables );
    startValues.resize( numberOfVariables );
    for( unsigned int whichIndex( 0 );
         whichIndex < numberOfVariables;
         ++whichIndex )
    {
      SetStartSystem( whichIndex,
                      lowerEndOfStartValues,
                      upperEndOfStartValues );
    }
  }

  // This evaluates the target system and places the values in
  // destinationVector.
  inline void
  PolynomialGradientTargetSystem::HomotopyContinuationSystemValues(
                   std::vector< std::complex< double > > solutionConfiguration,
                     std::vector< std::complex< double > >& destinationVector )
  {
    destinationVector.resize( targetSystem.size() );
    for( unsigned int whichIndex( 0 );
         whichIndex < targetSystem.size();
         ++whichIndex )
    {
      destinationVector[ whichIndex ]
      = targetSystem[ whichIndex ]( solutionConfiguration );
    }
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PolynomialGradientTargetSystem::HomotopyContinuationSystemValues("
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

  // This evaluates the derivatives of the target system and places the
  // values in destinationMatrix.
  inline void
  PolynomialGradientTargetSystem::HomotopyContinuationSystemGradients(
                   std::vector< std::complex< double > > solutionConfiguration,
      std::vector< std::vector< std::complex< double > > >& destinationMatrix )
  {
    destinationMatrix.resize( targetHessian.size() );
    for( unsigned int constraintIndex( 0 );
         constraintIndex < targetHessian.size();
         ++constraintIndex )
    {
      destinationMatrix[ constraintIndex ].resize(
                                     targetHessian[ constraintIndex ].size() );
      for( unsigned int variableIndex( 0 );
           variableIndex < targetHessian[ constraintIndex ].size();
           ++variableIndex )
      {
        destinationMatrix[ constraintIndex ][ variableIndex ]
        = targetHessian[ constraintIndex ][ variableIndex ](
                                                       solutionConfiguration );
      }
    }
  }

  // This fills targetSystem from potentialPolynomial.
  inline void PolynomialGradientTargetSystem::PreparePolynomialGradient(
                                     PolynomialSum const& potentialPolynomial )
  {
    targetSystem.assign( numberOfVariables,
                                     PolynomialSum() );
    std::vector< PolynomialTerm > const& polynomialForGradient(
                                       potentialPolynomial.PolynomialTerms() );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PolynomialGradientTargetSystem::PreparePolynomialGradient()"
    <<  " called. targetSystem.size() set to "
    << targetSystem.size()
    <<  ", polynomialForGradient.size() = " << polynomialForGradient.size();
    std::cout << std::endl;*/

    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfVariables;
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
          targetSystem[ fieldIndex ].PolynomialTerms().push_back(
                                  whichTerm->PartialDerivative( fieldIndex ) );
        }
      }
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "targetSystem =" << std::endl;
    for( std::vector< PolynomialSum >::iterator
         whichTadpole( targetSystem.begin() );
         whichTadpole < targetSystem.end();
         ++whichTadpole )
    {
      std::cout << std::endl << whichTadpole->AsDebuggingString();
    }
    std::cout << std::endl;
    std::cout << std::endl;*/
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALGRADIENTTARGETSYSTEM_HPP_ */
