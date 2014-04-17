/*
 * SumOfProductOfPolynomialSums.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SUMOFPRODUCTOFPOLYNOMIALSUMS_HPP_
#define SUMOFPRODUCTOFPOLYNOMIALSUMS_HPP_

#include "../CommonIncludes.hpp"
#include "ProductOfPolynomialSums.hpp"

namespace VevaciousPlusPlus
{

  class SumOfProductOfPolynomialSums
  {
  public:
    SumOfProductOfPolynomialSums();
    SumOfProductOfPolynomialSums(
                              SumOfProductOfPolynomialSums const& copySource );
    virtual
    ~SumOfProductOfPolynomialSums();


    // This returns the sum of operator() on each of the
    // ProductOfPolynomialSums in sumVector.
    std::complex< double > operator()(
       std::vector< std::complex< double > > const& fieldConfiguration ) const;

    // This pushes productOfPolynomialSums into sumVector.
    void AddToSum( ProductOfPolynomialSums const& productOfPolynomialSums )
    { sumVector.push_back( productOfPolynomialSums ); }

    // This is mainly for debugging.
    std::string AsString( std::vector< std::string > const& fieldNames ) const;


  protected:
    std::vector< ProductOfPolynomialSums > sumVector;
  };




  // This returns the sum of operator() on each of the
  // ProductOfPolynomialSums in sumVector.
  inline std::complex< double > SumOfProductOfPolynomialSums::operator()(
        std::vector< std::complex< double > > const& fieldConfiguration ) const
  {
    std::complex< double > returnValue( 0.0,
                                        0.0 );
    for( std::vector< ProductOfPolynomialSums >::const_iterator
         whichProduct( sumVector.begin() );
         whichProduct < sumVector.end();
         ++whichProduct )
    {
      returnValue += (*whichProduct)( fieldConfiguration );
    }
    return returnValue;
  }

  // This is mainly for debugging.
  inline std::string SumOfProductOfPolynomialSums::AsString(
                           std::vector< std::string > const& fieldNames ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << "(";
    for( std::vector< ProductOfPolynomialSums >::const_iterator
         whichProduct( sumVector.begin() );
         whichProduct < sumVector.end();
         ++whichProduct )
    {
      if( whichProduct != sumVector.begin() )
      {
        stringBuilder << "+";
      }
      stringBuilder << " " << whichProduct->AsString( fieldNames );
    }
    stringBuilder << " )";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* SUMOFPRODUCTOFPOLYNOMIALSUMS_HPP_ */
