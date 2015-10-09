/*
 * PolynomialSum.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALSUM_HPP_
#define POLYNOMIALSUM_HPP_

#include "CommonIncludes.hpp"
#include "PolynomialTerm.hpp"

namespace VevaciousPlusPlus
{
  class PolynomialSum
  {
  public:
    PolynomialSum();
    PolynomialSum( PolynomialSum const& copySource );
    virtual ~PolynomialSum();


    // This evaluates the sum of polynomial terms for the given field values,
    // with the funtionoids at their values from their last scale update.
    double operator()( std::vector< double > const& fieldConfiguration ) const;

    // This evaluates the sum of polynomial terms for the given field values,
    // with the funtionoids evaluated at the natural exponent of
    // logarithmOfScale.
    double operator()( std::vector< double > const& fieldConfiguration,
                       double const logarithmOfScale ) const;

    // This evaluates the sum of polynomial terms for the given field values.
    std::complex< double > operator()(
       std::vector< std::complex< double > > const& fieldConfiguration ) const;

    std::vector< PolynomialTerm > const& PolynomialTerms() const
    { return polynomialTerms; }

    std::vector< PolynomialTerm >& PolynomialTerms(){ return polynomialTerms; }

    // This returns the highest sum of field exponents of all the terms in
    // polynomialTerms.
    unsigned int HighestFieldPower() const;

    // This prints the polynomial with the current values of the coefficients
    // and the powers of the fields with names given by fieldNames to a string
    // and returns it.
    std::string AsStringAtCurrentScale(
                          std::vector< std::string > const& fieldNames ) const;

    // This is mainly for debugging:
    std::string AsDebuggingString() const;

    // This returns a string that should be valid Python assuming that the
    // field configuration is given as an array called "fv".
    std::string AsPython() const;


  protected:
    std::vector< PolynomialTerm > polynomialTerms;
  };




  // This evaluates the sum of polynomial terms for the given field values,
  // with the funtionoids at their values from their last scale update.
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

  // This evaluates the sum of polynomial terms for the given field values,
  // with the funtionoids evaluated at the natural exponent of
  // logarithmOfScale.
  inline double PolynomialSum::operator()(
                               std::vector< double > const& fieldConfiguration,
                                          double const logarithmOfScale ) const
  {
    double returnValue( 0.0 );
    for( std::vector< PolynomialTerm >::const_iterator
         whichTerm( polynomialTerms.begin() );
         whichTerm < polynomialTerms.end();
         ++whichTerm )
    {
      returnValue += (*whichTerm)( fieldConfiguration,
                                   logarithmOfScale );
    }
    return returnValue;
  }

  // This evaluates the sum of polynomial terms for the given field values.
  inline std::complex< double > PolynomialSum::operator()(
        std::vector< std::complex< double > > const& fieldConfiguration ) const
  {
    std::complex< double > returnValue( 0.0,
                                        0.0 );
    for( std::vector< PolynomialTerm >::const_iterator
         whichTerm( polynomialTerms.begin() );
         whichTerm < polynomialTerms.end();
         ++whichTerm )
    {
      returnValue += (*whichTerm)( fieldConfiguration );
    }
    return returnValue;
  }

  // This returns the highest sum of field exponents of all the terms in
  // polynomialTerms.
  inline unsigned int PolynomialSum::HighestFieldPower() const
  {
    unsigned int highestPower( 0 );
    for( std::vector< PolynomialTerm >::const_iterator
         whichTerm( polynomialTerms.begin() );
         whichTerm < polynomialTerms.end();
         ++whichTerm )
    {
      if( whichTerm->FieldPower() > highestPower )
      {
        highestPower = whichTerm->FieldPower();
      }
    }
    return highestPower;
  }

  // This prints the polynomial with the current values of the coefficients
  // and the powers of the fields with names given by fieldNames to a string
  // and returns it.
  inline std::string PolynomialSum::AsStringAtCurrentScale(
                        std::vector< std::string > const& fieldNames ) const
  {
    std::string returnString( "" );
    for( std::vector< PolynomialTerm >::const_iterator
         whichTerm( polynomialTerms.begin() );
         whichTerm < polynomialTerms.end();
         ++whichTerm )
    {
      returnString.append( whichTerm->AsStringAtCurrentScale( fieldNames,
                                      ( whichTerm != polynomialTerms.begin() ),
                                                              true ) );
    }
    return returnString;
  }

  // This is mainly for debugging:
  inline std::string PolynomialSum::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream << "PolynomialSum =" << std::endl;
    for( std::vector< PolynomialTerm >::const_iterator
         whichTerm( polynomialTerms.begin() );
         whichTerm < polynomialTerms.end();
         ++whichTerm )
    {
      returnStream << whichTerm->AsDebuggingString() << std::endl;
    }
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALSUM_HPP_ */
