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
                     std::vector< double > const& firstAndSecondDerivatives,
                     double const radialValue );


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
                        std::vector< double > const& firstAndSecondDerivatives,
                                          double const radialValue )
  {
    double const auxiliaryValue( auxiliaryAndFirstDerivative[ 0 ] );
    double const auxiliaryDerivative( auxiliaryAndFirstDerivative[ 1 ] );
    double firstDerivativeValue;
    double fieldDerivativeSquared( 0.0 );
    double fieldFirstDotSecondDerivatives( 0.0 );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      firstDerivativeValue = firstDerivatives( auxiliaryValue );
      fieldDerivativeSquared
      += ( firstDerivativeValue * firstDerivativeValue );
      fieldFirstDotSecondDerivatives
      += ( firstDerivativeValue * secondDerivatives( auxiliaryValue ) );
    }
    firstAndSecondDerivatives[ 0 ] = auxiliaryDerivative;
    firstAndSecondDerivatives[ 1 ]
     = ( ( ( potentialDerivative( auxiliaryValue )
             - ( fieldFirstDotSecondDerivatives
                 * auxiliaryDerivative * auxiliaryDerivative ) )
           / fieldDerivativeSquared )
         - ( ( dampingFactor * auxiliaryDerivative ) / radialValue ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* BUBBLEPROFILER_HPP_ */
