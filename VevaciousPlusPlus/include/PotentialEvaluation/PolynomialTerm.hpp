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

    // This returns true if the field with index fieldIndex has a non-zero
    // power.
    bool NonZeroDerivative( unsigned int const fieldIndex ) const;

    // This returns a PolynomialTerm that is the partial derivative with
    // respect to the field with index fieldIndex.
    PolynomialTerm PartialDerivative( unsigned int const fieldIndex ) const;

    // This removes any scale dependence from this PolynomialTerm, so its
    // coefficient will remain fixed unless it is multiplied by more
    // ParameterFunctionoid pointers.
    void RemoveScaleDependence(){ functionoidProduct.clear(); }

    // This is mainly for debugging.
    std::string AsDebuggingString() const;

    // This prints the polynomial with the current value of coefficientConstant
    // and the powers of the fields with names given by fieldNames to a string
    // and returns it.
    std::string AsStringAtCurrentScale(
                                  std::vector< std::string > const& fieldNames,
                              bool const prependPlusIfPositive = false ) const;


  protected:
    bool isValid;
    double coefficientConstant;
    std::vector< unsigned int > fieldProductByIndex;
    std::vector< unsigned int > fieldPowersByIndex;
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

  // This adds runningParameter to the set of functionoids which multiply
  // coefficientConstant to form the scale-dependent coefficient.
  inline void
  PolynomialTerm::MultiplyBy( ParameterFunctionoid* const runningParameter,
                              unsigned int const powerInt )
  {
    if( runningParameter == NULL )
    {
      isValid = false;
    }
    else
    {
      functionoidProduct.insert( functionoidProduct.end(),
                                 powerInt,
                                 runningParameter );
    }
  }

  // This raises the power of the field given by fieldIndex by the number
  // given by powerInt.
  inline void PolynomialTerm::RaiseFieldPower( unsigned int const fieldIndex,
                                               unsigned int const powerInt )
  {
    fieldProductByIndex.insert( fieldProductByIndex.end(),
                                powerInt,
                                fieldIndex );
    if( fieldPowersByIndex.size() <= fieldIndex )
    {
      fieldPowersByIndex.resize( ( fieldIndex + 1 ),
                                 0 );
    }
    fieldPowersByIndex[ fieldIndex ] += powerInt;
  }

  // This resets the PolynomialTerm to be as if freshly constructed.
  inline void PolynomialTerm::ResetValues()
  {
    isValid = true;
    coefficientConstant =  1.0;
    fieldProductByIndex.clear();
    fieldPowersByIndex.clear();
    functionoidProduct.clear();
  }

  // This returns true if the field with index fieldIndex has a non-zero
  // power.
  inline bool
  PolynomialTerm::NonZeroDerivative( unsigned int const fieldIndex ) const
  {
    return ( ( fieldPowersByIndex.size() > fieldIndex )
             &&
             ( fieldPowersByIndex[ fieldIndex ] > 0 ) );
  }

  // This returns a PolynomialTerm that is the partial derivative with
  // respect to the field with index fieldIndex.
  inline PolynomialTerm
  PolynomialTerm::PartialDerivative( unsigned int const fieldIndex ) const
  {
    if( !NonZeroDerivative( fieldIndex ) )
    {
      throw std::invalid_argument( "PolynomialTerm does not have non-zero"
                                  " derivative with respect to field index!" );
    }
    PolynomialTerm returnTerm( *this );
    returnTerm.coefficientConstant *= fieldPowersByIndex[ fieldIndex ];
    returnTerm.fieldPowersByIndex[ fieldIndex ] -= 1;
    returnTerm.fieldProductByIndex.clear();
    for( unsigned int whichField( 0 );
         whichField < fieldPowersByIndex.size();
         ++whichField )
    {
        returnTerm.fieldProductByIndex.insert(
                                          returnTerm.fieldProductByIndex.end(),
                                   returnTerm.fieldPowersByIndex[ whichField ],
                                               whichField );
    }
    return returnTerm;
  }

  // This prints the polynomial with the current value of coefficientConstant
  // and the powers of the fields with names given by fieldNames to a string
  // and returns it.
  inline std::string PolynomialTerm::AsStringAtCurrentScale(
                                  std::vector< std::string > const& fieldNames,
                                       bool const prependPlusIfPositive ) const
  {
    if( fieldNames.size() < fieldPowersByIndex.size() )
    {
      throw std::out_of_range( "Not enough field names for PolynomialTerm!" );
    }
    std::stringstream returnStream;
    if( prependPlusIfPositive
        &&
        ( coefficientConstant >= 0.0 ) )
    {
      returnStream << '+';
    }
    returnStream << coefficientConstant;
    for( unsigned int whichField( 0 );
         whichField < fieldPowersByIndex.size();
         ++whichField )
    {
      returnStream << " * " << fieldNames[ whichField ] << "^"
      << fieldPowersByIndex[ whichField ];
    }
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALTERM_HPP_ */
