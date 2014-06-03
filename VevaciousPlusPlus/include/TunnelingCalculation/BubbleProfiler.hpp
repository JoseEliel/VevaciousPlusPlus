/*
 * BubbleProfiler.hpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BUBBLEPROFILER_HPP_
#define BUBBLEPROFILER_HPP_

#include "../CommonIncludes.hpp"
#include "../PotentialEvaluation/SimplePolynomial.hpp"

namespace VevaciousPlusPlus
{

  class BubbleProfiler
  {
  public:
    BubbleProfiler( SimplePolynomial const& potentialApproximation,
                   std::vector< SimplePolynomial > const& fieldPathDerivatives,
                    unsigned int const dampingFactor );
    virtual
    ~BubbleProfiler();


    // This is in the form required for the Boost odeint package.
    void operator()( std::vector< double > const& auxiliaryAndFirstDerivative,
                     std::vector< double >& firstAndSecondDerivatives,
                     double const radialValue );

    // This sets initialConditions by a Euler step assuming that near r = 0,
    // p goes as p_0 + p_2 r^2 (as the bubble should have smooth fields at its
    // center); hence d^2p/dr^2 at r=0 is
    // ( dV/dp ) / ( ( 1 + 2 * dampingFactor ) |df/dp|^2 ).
    void DoFirstStep( std::vector< double >& initialConditions,
                      double const currentAuxiliary,
                      double const initialIntegrationStep );


  protected:
    SimplePolynomial potentialDerivative;
    unsigned int const numberOfFields;
    std::vector< SimplePolynomial > const& firstDerivatives;
    std::vector< SimplePolynomial > secondDerivatives;
    double const dampingFactor;
  };




  // This is in the form required for the Boost odeint package.
  inline void BubbleProfiler::operator()(
                      std::vector< double > const& auxiliaryAndFirstDerivative,
                              std::vector< double >& firstAndSecondDerivatives,
                                          double const radialValue )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "BubbleProfiler::operator()( { " << auxiliaryAndFirstDerivative[ 0 ]
    << ", " << auxiliaryAndFirstDerivative[ 1 ] << " }, { "
    << firstAndSecondDerivatives[ 0 ]
    << ", " << firstAndSecondDerivatives[ 1 ] << " }, " << radialValue
    << " ) called.";
    std::cout << std::endl;/**/
    double const auxiliaryValue( auxiliaryAndFirstDerivative[ 0 ] );
    double const auxiliaryDerivative( auxiliaryAndFirstDerivative[ 1 ] );
    double firstDerivativeValue( 0.0 );
    double fieldDerivativeSquared( 0.0 );
    double fieldFirstDotSecondDerivatives( 0.0 );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      firstDerivativeValue = firstDerivatives[ fieldIndex ]( auxiliaryValue );
      fieldDerivativeSquared
      += ( firstDerivativeValue * firstDerivativeValue );
      fieldFirstDotSecondDerivatives
      += ( firstDerivativeValue
           * secondDerivatives[ fieldIndex ]( auxiliaryValue ) );
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "1st derivative [" << fieldIndex << "] = " << firstDerivativeValue
      << std::endl
      << "2nd derivative [" << fieldIndex << "] = "
      << secondDerivatives[ fieldIndex ]( auxiliaryValue );
      std::cout << std::endl;/**/
    }
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "fieldDerivativeSquared = " << fieldDerivativeSquared
    << ", fieldFirstDotSecondDerivatives = " << fieldFirstDotSecondDerivatives;
    std::cout << std::endl;/**/
    firstAndSecondDerivatives[ 0 ] = auxiliaryDerivative;
    firstAndSecondDerivatives[ 1 ]
     = ( ( ( potentialDerivative( auxiliaryValue )
             - ( fieldFirstDotSecondDerivatives
                 * auxiliaryDerivative * auxiliaryDerivative ) )
           / fieldDerivativeSquared )
         - ( ( dampingFactor * auxiliaryDerivative ) / radialValue ) );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "firstAndSecondDerivatives = { " << firstAndSecondDerivatives[ 0 ]
    << ", " << firstAndSecondDerivatives[ 1 ] << " }";
    std::cout << std::endl;/**/
  }

  // This sets initialConditions by a Euler step assuming that near r = 0,
  // p goes as p_0 + p_2 r^2 (as the bubble should have smooth fields at its
  // center); hence d^2p/dr^2 at r = 0 is
  // ( dV/dp ) / ( ( 1 + 2 * dampingFactor ) |df/dp|^2 ).
  inline void
  BubbleProfiler::DoFirstStep( std::vector< double >& initialConditions,
                               double const currentAuxiliary,
                               double const initialIntegrationStep )
  {
    initialConditions[ 0 ] = currentAuxiliary;
    double firstDerivativeValue;
    double fieldDerivativeSquared( 0.0 );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      firstDerivativeValue
      = firstDerivatives[ fieldIndex ]( currentAuxiliary );
      fieldDerivativeSquared
      += ( firstDerivativeValue * firstDerivativeValue );
    }
    initialConditions[ 1 ] = ( ( potentialDerivative( currentAuxiliary )
                                 * initialIntegrationStep )
                               / ( ( 1.0 + dampingFactor + dampingFactor )
                                   * fieldDerivativeSquared ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* BUBBLEPROFILER_HPP_ */
