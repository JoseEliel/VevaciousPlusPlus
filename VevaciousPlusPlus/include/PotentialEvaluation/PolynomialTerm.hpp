/*
 * PolynomialTerm.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALTERM_HPP_
#define POLYNOMIALTERM_HPP_

#include "../StandardIncludes.hpp"
#include "ParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialTerm
  {
  public:
    PolynomialTerm( bool const coefficientIsPositive );
    virtual
    ~PolynomialTerm();


    // This returns a flag that can be false if the coefficient got multiplied
    // by a NULL functionoid.
    bool IsValid() const { return isValid; }

    // This multiplies coefficientConstant by multiplicationFactor.
    void MultiplyBy( double const multiplicationFactor );

    // This adds runningParameter to the set of functionoids which multiply
    // coefficientConstant to form the scale-dependent coefficient.
    void MultiplyBy( ParameterFunctionoid* const runningParameter,
                     int const powerInt );



  protected:
    bool isValid;
    double coefficientConstant;
  };

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALTERM_HPP_ */
