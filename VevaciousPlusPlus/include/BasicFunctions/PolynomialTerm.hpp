/*
 * PolynomialTerm.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALTERM_HPP_
#define POLYNOMIALTERM_HPP_

#include "CommonIncludes.hpp"
#include "PotentialEvaluation/ParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialTerm : public BOL::BasicObserver
  {
  public:
    PolynomialTerm();
    PolynomialTerm( PolynomialTerm const& copySource );
    virtual ~PolynomialTerm();


    // This multiplies the relevant field values with the coefficient and the
    // values from the functionoids at the scale of the last call of
    // UpdateForUpdatedFunctionoids.
    double operator()( std::vector< double > const& fieldConfiguration ) const
    { return FieldProduct( currentScaleTotalCoefficient,
                           fieldConfiguration ); }

    // This multiplies the relevant field values with the coefficient and the
    // values from the functionoids at the scale given by the natural exponent
    // of logarithmOfScale.
    double operator()( std::vector< double > const& fieldConfiguration,
                       double const logarithmOfScale ) const
    { return FunctionoidProduct( FieldProduct( coefficientConstant,
                                               fieldConfiguration ),
                                 logarithmOfScale ); }

    // This multiplies the relevant field values with the coefficient and the
    // values from the functionoids at the scale of the last call of
    // UpdateForUpdatedFunctionoids.
    std::complex< double > operator()(
       std::vector< std::complex< double > > const& fieldConfiguration ) const;

    // This updates currentScaleTotalCoefficient to be coefficientConstant
    // multiplied by the product of the operator() values from
    // functionoidProduct.
    void respondToObservedSignal(){ currentScaleTotalCoefficient
                                 = FunctionoidProduct( coefficientConstant ); }

    // This returns a flag that can be false if the coefficient got multiplied
    // by a NULL functionoid.
    bool IsValid() const { return isValid; }

    // This raises the power of the field given by fieldIndex by the number
    // given by powerInt.
    void RaiseFieldPower( size_t const fieldIndex,
                          size_t const powerInt );

    // This multiplies coefficientConstant by multiplicationFactor.
    void MultiplyBy( double const multiplicationFactor );

    // This adds runningParameter to the set of functionoids which multiply
    // coefficientConstant to form the scale-dependent coefficient.
    void MultiplyBy( ParameterFunctionoid* const runningParameter,
                     size_t const powerInt );

    // This resets the PolynomialTerm to be as if freshly constructed.
    void ResetValues();

    // This returns true if the field with index fieldIndex has a non-zero
    // power.
    bool NonZeroDerivative( size_t const fieldIndex ) const
    { return ( ( fieldPowersByIndex.size() > fieldIndex )
               &&
               ( fieldPowersByIndex[ fieldIndex ] > 0 ) ); }

    // This returns a PolynomialTerm that is the partial derivative with
    // respect to the field with index fieldIndex.
    PolynomialTerm PartialDerivative( size_t const fieldIndex ) const;

    // This returns the sum of the powers of the fields.
    size_t FieldPower() const{ return fieldProductByIndex.size(); }

    // This removes any scale dependence from this PolynomialTerm, so its
    // coefficient will remain fixed unless it is multiplied by more
    // ParameterFunctionoid pointers.
    void RemoveScaleDependence(){ functionoidProduct.clear(); }

    // This prints the polynomial with the current value of coefficientConstant
    // and the powers of the fields with names given by fieldNames to a string
    // and returns it.
    std::string AsStringAtCurrentScale(
                                  std::vector< std::string > const& fieldNames,
                                        bool const followsOtherTerm = false,
                                        bool const emptyIfZero = false ) const;

    // This is mainly for debugging.
    std::string AsDebuggingString() const;

    // This returns a string that should be valid Python assuming that the
    // field configuration is given as an array called "fv".
    std::string AsPython() const;


  protected:
    bool isValid;
    double coefficientConstant;
    std::vector< size_t > fieldProductByIndex;
    std::vector< size_t > fieldPowersByIndex;
    std::vector< ParameterFunctionoid* > functionoidProduct;
    double currentScaleTotalCoefficient;

    // This returns doubleToMultiply multiplied by the field product using the
    // values in fieldConfiguration.
    double const FieldProduct( double doubleToMultiply,
                       std::vector< double > const& fieldConfiguration ) const;

    // This returns doubleToMultiply multiplied by the funtionoids in
    // functionoidProduct using the values from the last scale update.
    double const FunctionoidProduct( double doubleToMultiply ) const;

    // This returns doubleToMultiply multiplied by the funtionoids in
    // functionoidProduct at the scale given by the natural exponent of
    // logarithmOfScale.
    double const FunctionoidProduct( double doubleToMultiply,
                                     double const logarithmOfScale ) const;
  };




  // This multiplies the relevant field values with the coefficient and the
  // values from the functionoids at the scale of the last call of
  // UpdateForUpdatedFunctionoids.
  inline std::complex< double > PolynomialTerm::operator()(
        std::vector< std::complex< double > > const& fieldConfiguration ) const
  {
    std::complex< double > returnValue( currentScaleTotalCoefficient,
                                        0.0 );
    for( std::vector< size_t >::const_iterator
         whichField( fieldProductByIndex.begin() );
         whichField < fieldProductByIndex.end();
         ++whichField )
    {
      returnValue *= fieldConfiguration[ *whichField ];
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
                              size_t const powerInt )
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
      runningParameter->registerObserver( this );
    }
  }

  // This raises the power of the field given by fieldIndex by the number
  // given by powerInt.
  inline void PolynomialTerm::RaiseFieldPower( size_t const fieldIndex,
                                               size_t const powerInt )
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
    coefficientConstant = 1.0;
    fieldProductByIndex.clear();
    fieldPowersByIndex.clear();
    functionoidProduct.clear();
    currentScaleTotalCoefficient = 1.0;
  }

  // This returns a PolynomialTerm that is the partial derivative with
  // respect to the field with index fieldIndex.
  inline PolynomialTerm
  PolynomialTerm::PartialDerivative( size_t const fieldIndex ) const
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
    for( size_t whichField( 0 );
         whichField < fieldPowersByIndex.size();
         ++whichField )
    {
        returnTerm.fieldProductByIndex.insert(
                                          returnTerm.fieldProductByIndex.end(),
                                   returnTerm.fieldPowersByIndex[ whichField ],
                                               whichField );
    }
    returnTerm.respondToObservedSignal();
    return returnTerm;
  }

  // This returns doubleToMultiply multiplied by the field product using the
  // values in fieldConfiguration.
  inline double const PolynomialTerm::FieldProduct( double doubleToMultiply,
                        std::vector< double > const& fieldConfiguration ) const
  {
    for( std::vector< size_t >::const_iterator
         whichField( fieldProductByIndex.begin() );
         whichField < fieldProductByIndex.end();
         ++whichField )
    {
      doubleToMultiply *= fieldConfiguration[ *whichField ];
    }
    return doubleToMultiply;
  }

  // This returns doubleToMultiply multiplied by the funtionoids in
  // functionoidProduct using the values from the last scale update.
  inline double const
  PolynomialTerm::FunctionoidProduct( double doubleToMultiply ) const
  {
    for( std::vector< ParameterFunctionoid* >::const_iterator
         whichFunctionoid( functionoidProduct.begin() );
         whichFunctionoid < functionoidProduct.end();
         ++whichFunctionoid )
    {
      doubleToMultiply *= (*(*whichFunctionoid))();
    }
    return doubleToMultiply;
  }

  // This returns doubleToMultiply multiplied by the funtionoids in
  // functionoidProduct at the scale given by the natural exponent of
  // logarithmOfScale.
  inline double const
  PolynomialTerm::FunctionoidProduct( double doubleToMultiply,
                                      double const logarithmOfScale ) const
  {
    for( std::vector< ParameterFunctionoid* >::const_iterator
         whichFunctionoid( functionoidProduct.begin() );
         whichFunctionoid < functionoidProduct.end();
         ++whichFunctionoid )
    {
      doubleToMultiply *= (*(*whichFunctionoid))( logarithmOfScale );
    }
    return doubleToMultiply;
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALTERM_HPP_ */
