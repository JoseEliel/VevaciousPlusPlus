/*
 * PolynomialGradientTargetSystem.hpp
 *
 *  Created on: Apr 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALGRADIENTTARGETSYSTEM_HPP_
#define POLYNOMIALGRADIENTTARGETSYSTEM_HPP_

#include "CommonIncludes.hpp"
#include "Eigen/Dense"
#include "HomotopyContinuationTargetSystem.hpp"
#include "SlhaManagement/SlhaUpdatePropagator.hpp"
#include "BasicFunctions/PolynomialTerm.hpp"
#include "BasicFunctions/PolynomialSum.hpp"
#include "BasicFunctions/ProductOfPolynomialSums.hpp"
#include "BasicFunctions/SumOfProductOfPolynomialSums.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialGradientTargetSystem :
                                        public HomotopyContinuationTargetSystem
  {
  public:
    PolynomialGradientTargetSystem( PolynomialSum const& potentialPolynomial,
                                    size_t const numberOfVariables,
                                    SlhaUpdatePropagator& previousPropagator,
                            std::vector< size_t > const& fieldsAssumedPositive,
                            std::vector< size_t > const& fieldsAssumedNegative,
                                bool const treeLevelMinimaOnlyAsValidSolutions,
                             double const assumedPositiveOrNegativeTolerance );
    virtual ~PolynomialGradientTargetSystem();


    // This fills targetSystem, startSystem, targetHessian, and startHessian
    // appropriately. The starting system is set with
    // SetStartSystem( variableIndex,
    //                 slhaManager.LowestBlockScale(),
    //                 slhaManager.HighestBlockScale() )
    // for variableIndex from 0 to ( numberOfVariables - 1 ). SetStartSystem
    // is public, so special ranges for particular variables can be given to
    // over-write the values given by this function.
    virtual void UpdateSelfForNewSlha( SlhaManager const& slhaManager );

    std::vector< PolynomialSum > const& TargetSystem() const
    { return targetSystem; }

    std::vector< PolynomialSum >& TargetSystem(){ return targetSystem; }

    std::vector< std::vector< PolynomialSum > > const& TargetHessian() const
    { return targetHessian; }

    std::vector< ProductOfPolynomialSums > const& StartSystem() const
    { return startSystem; }

    std::vector< std::vector< std::complex< double > > > const&
    StartValues() const{ return startValues; }

    std::vector< std::vector< SumOfProductOfPolynomialSums > > const&
    StartHessian() const{ return startHessian; }

    // This evaluates the target system and places the values in
    // destinationVector.
    virtual void HomotopyContinuationSystemValues(
                   std::vector< std::complex< double > > solutionConfiguration,
              std::vector< std::complex< double > >& destinationVector ) const;

    // This is the real-valued version of the above.
    virtual void HomotopyContinuationSystemValues(
                                   std::vector< double > solutionConfiguration,
                              std::vector< double >& destinationVector ) const;

    // This evaluates the derivatives of the target system and places the
    // values in destinationMatrix.
    virtual void HomotopyContinuationSystemGradients(
                   std::vector< std::complex< double > > solutionConfiguration,
        std::vector< std::vector< std::complex< double > > >& destinationMatrix
                                                                       ) const;

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
    void SetStartSystem( size_t const variableIndex,
                         double const lowerEndOfStartValues,
                         double const upperEndOfStartValues );


  protected:
    PolynomialSum const& potentialPolynomial;
    size_t const numberOfVariables;
    size_t numberOfFields;
    std::vector< PolynomialSum > targetSystem;
    size_t highestDegreeOfTargetSystem;
    double minimumScale;
    double maximumScale;
    std::vector< ProductOfPolynomialSums > startSystem;
    std::vector< std::vector< std::complex< double > > > startValues;
    std::vector< std::vector< std::complex< double > > > validSolutions;
    std::vector< std::vector< PolynomialSum > > targetHessian;
    std::vector< std::vector< SumOfProductOfPolynomialSums > > startHessian;
    std::vector< size_t > const& fieldsAssumedPositive;
    std::vector< size_t > const& fieldsAssumedNegative;
    bool treeLevelMinimaOnlyAsValidSolutions;
    double const assumedPositiveOrNegativeTolerance;


    // This vetoes a homotopy continuation solution if any of the fields with
    // index in fieldsAssumedPositive are negative (allowing for a small amount
    // of numerical jitter given by equationTolerance) or if any of the fields
    // with index in fieldsAssumedNegitive are positive (also allowing for a
    // small amount of numerical jitter), or if
    // treeLevelMinimaOnlyAsValidSolutions is true and the solution does not
    // correspond to a minimum (rather than just an extremum) of
    // potentialPolynomial.
    virtual bool AllowedSolution(
                    std::vector< double > const& solutionConfiguration ) const;

    // This fills targetSystem from potentialPolynomial.
    virtual void
    PreparePolynomialGradient( PolynomialSum const& potentialPolynomial );

    // This fills targetHessian from targetSystem.
    virtual void PreparePolynomialHessian();

    // This fills startHessian[ variableIndex ] from
    // startValues[ variableIndex ].
    void SetStartHessian( size_t const variableIndex );

    // This returns the 2nd derivative matrix of potentialPolynomial with
    // respect to its fields.
    virtual Eigen::MatrixXd FieldHessianOfPotential(
                    std::vector< double > const& solutionConfiguration ) const;
  };




  // This fills targetSystem, startSystem, targetHessian, and startHessian
  // appropriately. The starting system is set with
  // SetStartSystem( variableIndex,
  //                 slhaManager.LowestBlockScale(),
  //                 slhaManager.HighestBlockScale() )
  // for variableIndex from 0 to ( numberOfVariables - 1 ).
  inline void PolynomialGradientTargetSystem::UpdateSelfForNewSlha(
                                               SlhaManager const& slhaManager )
  {
    highestDegreeOfTargetSystem
    = ( potentialPolynomial.HighestFieldPower() - 1 );
    PreparePolynomialGradient( potentialPolynomial );
    PreparePolynomialHessian();
    startSystem.resize( numberOfVariables );
    startValues.resize( numberOfVariables );
    startHessian.resize( numberOfVariables );
    minimumScale = slhaManager.LowestBlockScale();
    maximumScale = slhaManager.HighestBlockScale();
    for( size_t whichIndex( 0 );
         whichIndex < numberOfVariables;
         ++whichIndex )
    {
      SetStartSystem( whichIndex,
                      minimumScale,
                      minimumScale );
    }
  }

  // This evaluates the target system and places the values in
  // destinationVector.
  inline void
  PolynomialGradientTargetSystem::HomotopyContinuationSystemValues(
                   std::vector< std::complex< double > > solutionConfiguration,
               std::vector< std::complex< double > >& destinationVector ) const
  {
    destinationVector.resize( targetSystem.size() );
    for( size_t whichIndex( 0 );
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
    for( size_t whichIndex( 0 );
         whichIndex < solutionConfiguration.size();
         ++whichIndex )
    {
      std::cout << " ( " << solutionConfiguration[ whichIndex ].real()
       << ", " << solutionConfiguration[ whichIndex ].imag() << " )";
    }
    std::cout << " }, ... ) set destinationVector to be { ";
    for( size_t whichIndex( 0 );
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

  // This evaluates the target system and places the values in
  // destinationVector.
  inline void
  PolynomialGradientTargetSystem::HomotopyContinuationSystemValues(
                                   std::vector< double > solutionConfiguration,
                               std::vector< double >& destinationVector ) const
  {
    destinationVector.resize( targetSystem.size() );
    for( size_t whichIndex( 0 );
         whichIndex < targetSystem.size();
         ++whichIndex )
    {
      destinationVector[ whichIndex ]
      = targetSystem[ whichIndex ]( solutionConfiguration );
    }
  }

  // This evaluates the derivatives of the target system and places the
  // values in destinationMatrix.
  inline void
  PolynomialGradientTargetSystem::HomotopyContinuationSystemGradients(
                   std::vector< std::complex< double > > solutionConfiguration,
        std::vector< std::vector< std::complex< double > > >& destinationMatrix
                                                                        ) const
  {
    destinationMatrix.resize( targetHessian.size() );
    for( size_t constraintIndex( 0 );
         constraintIndex < targetHessian.size();
         ++constraintIndex )
    {
      destinationMatrix[ constraintIndex ].resize(
                                     targetHessian[ constraintIndex ].size() );
      for( size_t variableIndex( 0 );
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

    for( size_t fieldIndex( 0 );
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
    << "potentialPolynomial =" << std::endl
    << potentialPolynomial.AsDebuggingString() << std::endl
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


  // This returns the 2nd derivative matrix of potentialPolynomial with
  // respect to its fields.
  inline Eigen::MatrixXd
  PolynomialGradientTargetSystem::FieldHessianOfPotential(
                     std::vector< double > const& solutionConfiguration ) const
  {
    Eigen::MatrixXd hessianMatrix( numberOfFields,
                                   numberOfFields );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfFields;
         ++rowIndex )
    {
      for( size_t columnIndex( 0 );
           columnIndex < numberOfFields;
           ++columnIndex )
      {
        hessianMatrix( rowIndex,
                       columnIndex )
        = targetHessian[ rowIndex ][ columnIndex ]( solutionConfiguration );
      }
    }
    return hessianMatrix;
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALGRADIENTTARGETSYSTEM_HPP_ */
