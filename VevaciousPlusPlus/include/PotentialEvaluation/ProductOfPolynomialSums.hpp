/*
 * ProductOfPolynomialSums.hpp
 *
 *  Created on: Apr 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PRODUCTOFPOLYNOMIALSUMS_HPP_
#define PRODUCTOFPOLYNOMIALSUMS_HPP_

#include "../StandardIncludes.hpp"
#include "PolynomialSum.hpp"

namespace VevaciousPlusPlus
{

  class ProductOfPolynomialSums
  {
  public:
    ProductOfPolynomialSums();
    virtual
    ~ProductOfPolynomialSums();


    // This returns the product of operator() on each of the PolynomialSums in
    // PolynomialSum.
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


  protected:
    std::vector< PolynomialSum > productVector;
  };




  // This returns the product of operator() on each of the PolynomialSums in
  // PolynomialSum.
  inline std::complex< double > ProductOfPolynomialSums::operator()(
        std::vector< std::complex< double > > const& fieldConfiguration ) const
  {
    double returnValue( 1.0 );
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

} /* namespace VevaciousPlusPlus */
#endif /* PRODUCTOFPOLYNOMIALSUMS_HPP_ */
