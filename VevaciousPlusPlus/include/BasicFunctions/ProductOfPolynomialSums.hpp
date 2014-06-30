/*
 * ProductOfPolynomialSums.hpp
 *
 *  Created on: Apr 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PRODUCTOFPOLYNOMIALSUMS_HPP_
#define PRODUCTOFPOLYNOMIALSUMS_HPP_

#include "CommonIncludes.hpp"
#include "PolynomialSum.hpp"

namespace VevaciousPlusPlus
{

  class ProductOfPolynomialSums
  {
  public:
    ProductOfPolynomialSums();
    ProductOfPolynomialSums( ProductOfPolynomialSums const& copySource );
    virtual
    ~ProductOfPolynomialSums();


    // This returns the product of operator() on each of the PolynomialSums in
    // productVector.
    std::complex< double > operator()(
       std::vector< std::complex< double > > const& fieldConfiguration ) const;

    // This pushes polynomialSum into productVector.
    void MultiplyBy( PolynomialSum const& polynomialSum )
    { productVector.push_back( polynomialSum ); }

    // This constructs f^fieldPower - startValue^fieldPower, where f is the
    // field with index fieldIndex, and pushes it into productVector.
    void MultiplyBy( unsigned int const fieldIndex,
                     unsigned int const fieldPower,
                     double const startValue );

    // This is mainly for debugging.
    std::string AsString( std::vector< std::string > const& fieldNames ) const;


  protected:
    std::vector< PolynomialSum > productVector;
  };




  // This returns the product of operator() on each of the PolynomialSums in
  // productVector.
  inline std::complex< double > ProductOfPolynomialSums::operator()(
        std::vector< std::complex< double > > const& fieldConfiguration ) const
  {
    std::complex< double > returnValue( 1.0,
                                        1.0 );
    for( std::vector< PolynomialSum >::const_iterator
         whichProduct( productVector.begin() );
         whichProduct < productVector.end();
         ++whichProduct )
    {
      returnValue *= (*whichProduct)( fieldConfiguration );
    }
    return returnValue;
  }

  // This constructs f^fieldPower - startValue^fieldPower, where f is the
  // field with index fieldIndex, and pushes it into productVector.
  inline void
  ProductOfPolynomialSums::MultiplyBy( unsigned int const fieldIndex,
                                       unsigned int const fieldPower,
                                       double const startValue )
  {
    PolynomialTerm fieldTerm;
    fieldTerm.RaiseFieldPower( fieldIndex,
                               fieldPower );
    PolynomialTerm constantTerm;
    constantTerm.MultiplyBy( -1.0 * pow( startValue,
                                         fieldPower ) );
    PolynomialSum constructedSum;
    constructedSum.PolynomialTerms().push_back( fieldTerm );
    constructedSum.PolynomialTerms().push_back( constantTerm );
    productVector.push_back( constructedSum );
  }

  // This is mainly for debugging.
  inline std::string ProductOfPolynomialSums::AsString(
                           std::vector< std::string > const& fieldNames ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << "[";
    for( std::vector< PolynomialSum >::const_iterator
         whichProduct( productVector.begin() );
         whichProduct < productVector.end();
         ++whichProduct )
    {
      if( whichProduct != productVector.begin() )
      {
        stringBuilder << " *";
      }
      stringBuilder << " ("
      << whichProduct->AsStringAtCurrentScale( fieldNames ) << ")";
    }
    stringBuilder << " ]";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* PRODUCTOFPOLYNOMIALSUMS_HPP_ */
