/*
 * OdeintBubbleDerivatives.hpp
 *
 *  Created on: Jun 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef ODEINTBUBBLEDERIVATIVES_HPP_
#define ODEINTBUBBLEDERIVATIVES_HPP_

#include "../CommonIncludes.hpp"
#include "../PotentialEvaluation/SimplePolynomial.hpp"
#include "BubbleRadialValueDescription.hpp"

namespace VevaciousPlusPlus
{

  class OdeintBubbleDerivatives
  {
  public:
    OdeintBubbleDerivatives(
                        PathFieldsAndPotential const& pathFieldsAndPotential );
    virtual
    ~OdeintBubbleDerivatives();

    // This puts the first and second derivatives based on
    // auxiliaryAndFirstDerivative into firstAndSecondDerivatives, in the form
    // required for the Boost odeint package.
    void operator()( std::vector< double > const& auxiliaryAndFirstDerivative,
                     std::vector< double >& firstAndSecondDerivatives,
                     double const radialValue );

    // This returns the first derivative of the polynomial approximation of the
    // potential function along the path with respect to the path auxiliary.
    double PotentialDerivative( double const auxiliaryValue ) const
    { return potentialDerivative( auxiliaryValue ); }


  protected:
    SimplePolynomial potentialDerivative;
    std::vector< SimplePolynomial > const& firstDerivatives;
    size_t const numberOfFields;
    std::vector< SimplePolynomial > secondDerivatives;
    double const dampingFactor;
  };




  // This is in the form required for the Boost odeint package.
  inline void OdeintBubbleDerivatives::operator()(
                      std::vector< double > const& auxiliaryAndFirstDerivative,
                              std::vector< double >& firstAndSecondDerivatives,
                                          double const radialValue )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "OdeintBubbleDerivatives::operator()( { "
    << auxiliaryAndFirstDerivative[ 0 ]
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

} /* namespace VevaciousPlusPlus */
#endif /* ODEINTBUBBLEDERIVATIVES_HPP_ */
