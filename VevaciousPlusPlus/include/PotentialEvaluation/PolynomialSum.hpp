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

    // This is mainly for debugging:
    std::string AsString() const;

  protected:
    std::vector< PolynomialTerm > polynomialTerms;
  };




  // This evaluates the sum of polynomial terms for the given field values.
  inline double PolynomialSum::operator()(
                        std::vector< double > const& fieldConfiguration ) const
  {
    double returnValue( 0.0 );
    for( std::vector< PolynomialTerm >::const_iterator
         whichTerm( polynomialTerms.begin() );
         whichTerm < polynomialTerms.end();
         ++whichTerm )
    {
      returnValue += (*whichTerm)( fieldConfiguration );
    }
    return returnValue;
  }

  // This is mainly for debugging:
  inline std::string PolynomialSum::AsString() const
  {
    std::stringstream returnStream;
    returnStream << "PolynomialSum =" << std::endl;
    for( std::vector< PolynomialTerm >::const_iterator
         whichTerm( polynomialTerms.begin() );
         whichTerm < polynomialTerms.end();
         ++whichTerm )
    {
      returnStream << whichTerm->AsString() << std::endl;
    }
    return std::string( returnStream.str() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALSUM_HPP_ */
