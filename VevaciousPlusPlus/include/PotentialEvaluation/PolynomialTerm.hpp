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
    PolynomialTerm();
    PolynomialTerm( PolynomialTerm const& copySource );
    virtual
    ~PolynomialTerm();


    // This multiplies the relevant field values with the coefficient and the
    // values from the functionoids.
    double operator()( std::vector< double > const& fieldConfiguration ) const;

    // This returns a flag that can be false if the coefficient got multiplied
    // by a NULL functionoid.
    bool IsValid() const { return isValid; }

    // This raises the power of the field given by fieldIndex by the number
    // given by powerInt.
    void RaiseFieldPower( unsigned int const fieldIndex,
                          unsigned int const powerInt );

    // This multiplies coefficientConstant by multiplicationFactor.
    void MultiplyBy( double const multiplicationFactor );

    // This adds runningParameter to the set of functionoids which multiply
    // coefficientConstant to form the scale-dependent coefficient.
    void MultiplyBy( ParameterFunctionoid* const runningParameter,
                     unsigned int const powerInt );

    // This resets the PolynomialTerm to be as if freshly constructed.
    void ResetValues();

    // This is mainly for debugging.
    std::string AsString() const;


  protected:
    bool isValid;
    double coefficientConstant;
    std::vector< unsigned int > fieldProductByIndex;
    std::vector< ParameterFunctionoid* > functionoidProduct;
  };




  // This multiplies the relevant field values with the coefficient and the
  // values from the functionoids.
  inline double PolynomialTerm::operator()(
                        std::vector< double > const& fieldConfiguration ) const
  {
    double returnValue( coefficientConstant );
    for( std::vector< unsigned int >::const_iterator
         whichField( fieldProductByIndex.begin() );
         whichField < fieldProductByIndex.end();
         ++whichField )
    {
      returnValue *= fieldConfiguration[ *whichField ];
    }
    for( std::vector< ParameterFunctionoid* >::const_iterator
         whichFunctionoid( functionoidProduct.begin() );
         whichFunctionoid < functionoidProduct.end();
         ++whichFunctionoid )
    {
      returnValue *= (*(*whichFunctionoid))();
    }
    return returnValue;
  }

  // This multiplies coefficientConstant by multiplicationFactor.
  inline void PolynomialTerm::MultiplyBy( double const multiplicationFactor )
  {
    coefficientConstant *= multiplicationFactor;
  }

  // This raises the power of the field given by fieldIndex by the number
  // given by powerInt.
  inline void PolynomialTerm::RaiseFieldPower( unsigned int const fieldIndex,
                                               unsigned int const powerInt )
  {
    for( unsigned int powerCount( 0 );
         powerCount < powerInt;
         ++powerCount )
    {
      fieldProductByIndex.push_back( fieldIndex );
    }
  }

  // This resets the PolynomialTerm to be as if freshly constructed.
  inline void PolynomialTerm::ResetValues()
  {
    isValid = true;
    coefficientConstant =  1.0;
    fieldProductByIndex.clear();
    functionoidProduct.clear();
  }

  // This is mainly for debugging.
  inline std::string PolynomialTerm::AsString() const
  {
    std::stringstream returnStream;
    returnStream
    << "isValid = " << isValid << std::endl
    << "coefficientConstant = " << coefficientConstant << std::endl
    << "fieldProductByIndex = {";
    for( std::vector< unsigned int >::const_iterator
         fieldIndex( fieldProductByIndex.begin() );
         fieldIndex < fieldProductByIndex.end();
         ++fieldIndex )
    {
      returnStream << " " << *fieldIndex;
    }
    returnStream
    << " }" << std::endl
    << "functionoidProduct = {";
    for( std::vector< ParameterFunctionoid* >::const_iterator
         functionoidPointer( functionoidProduct.begin() );
         functionoidPointer < functionoidProduct.end();
         ++functionoidPointer )
    {
      returnStream << " " << *functionoidPointer;
    }
    returnStream
    << " }" << std::endl;
    return std::string( returnStream.str() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALTERM_HPP_ */
