/*
 * PolynomialSum.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALSUM_HPP_
#define POLYNOMIALSUM_HPP_

#include  "../StandardIncludes.hpp"
#include "PolynomialTerm.hpp"

namespace VevaciousPlusPlus
{
  class PolynomialSum
  {
  public:
    PolynomialSum();
    virtual
    ~PolynomialSum();


    // This evaluates the sum of polynomial terms for the given field values.
    double operator()( std::vector< double > const& fieldConfiguration ) const;

    std::vector< PolynomialTerm >& PolynomialTerms(){ return polynomialTerms; }


  protected:
    std::vector< PolynomialTerm > polynomialTerms;
  };




  // This evaluates the sum of polynomial terms for the given field values.
  inline double PolynomialSum::operator()(
                        std::vector< double > const& fieldConfiguration ) const
  {
    double returnValue( 0.0 );
    for( std::vector< PolynomialTerm >::iterator
         whichTerm( polynomialTerms.begin() );
         whichTerm < polynomialTerms.end();
         ++whichTerm )
    {
      returnValue += (*whichTerm)( fieldConfiguration );
    }
    return returnValue;
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALSUM_HPP_ */
