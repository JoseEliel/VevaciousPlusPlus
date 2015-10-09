/*
 * ParametersAndFieldsProduct.hpp
 *
 *  Created on: Oct 8, 2015
 *      Author: bol
 */

#ifndef PARAMETERSANDFIELDSPRODUCT_HPP_
#define PARAMETERSANDFIELDSPRODUCT_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{
  // This class serves the purpose of turning a vector of field values into a
  // value of a Lagrangian term. Either it uses values for the Lagrangian
  // parameters stored from a previous lookup of the parameters for a
  // fixed-scale calculation, or it multiplies out selected values from a given
  // vector of Lagrangian parameters for the relevant scale.
  class ParametersAndFieldsProduct
  {
  public:
    ParametersAndFieldsProduct();
    ParametersAndFieldsProduct( ParametersAndFieldsProduct const& copySource );
    virtual ~ParametersAndFieldsProduct();


    // This sets the coefficient which multiplies the field values once for a
    // new parameter point, to save the product of parameters being
    // re-calculated many times for the same parameter point.
    void UpdateForFixedScale( std::vector< double > const& parameterValues )
    { totalCoefficientForFixedScale = ElementProduct( coefficientConstant,
                                                      parameterValues,
                                                      parameterIndices ); }

    // This multiplies the relevant field values with the coefficient and the
    // values from the Lagrangian parameters found in parameterValues.
    double operator()( std::vector< double > const& parameterValues,
                       std::vector< double > const& fieldConfiguration ) const
    { return ElementProduct( ElementProduct( coefficientConstant,
                                             fieldConfiguration,
                                             fieldProductByIndex ),
                             parameterValues,
                             parameterIndices ); }

    // This multiplies the relevant field values with the coefficient and the
    // values of the Lagrangian parameters from the last call of
    // UpdateForFixedScale.
    double operator()( std::vector< double > const& fieldConfiguration ) const
    { return ElementProduct( totalCoefficientForFixedScale,
                             fieldConfiguration,
                             fieldProductByIndex ); }

    // This raises the power of the field given by fieldIndex by the number
    // given by powerInt.
    void RaiseFieldPower( size_t const fieldIndex,
                          size_t const powerInt );

    // This multiplies coefficientConstant by multiplicationFactor.
    void MultiplyByConstant( double const multiplicationFactor )
    { coefficientConstant *= multiplicationFactor; }

    // This adds parameterIndex to the set of indices used to select the values
    // of Lagrangian parameters when forming the scale-dependent coefficient
    // given a vector of Lagrangian parameters evaluated at the relevant scale.
    void MultiplyByParameter( size_t const parameterIndex,
                              size_t const powerInt )
    { parameterIndices.insert( parameterIndices.end(),
                               powerInt,
                               parameterIndex ); }

    // This resets the ParametersAndFieldsProduct to be as if freshly
    // constructed.
    void ResetValues();

    // This returns true if the field with index fieldIndex has a non-zero
    // power.
    bool NonZeroDerivative( size_t const fieldIndex ) const
    { return ( ( fieldPowersByIndex.size() > fieldIndex )
               &&
               ( fieldPowersByIndex[ fieldIndex ] > 0 ) ); }

    // This returns a ParametersAndFieldsProduct that is the partial derivative
    // with respect to the field with index fieldIndex.
    ParametersAndFieldsProduct
    PartialDerivative( size_t const fieldIndex ) const;

    // This returns the sum of the powers of the fields.
    size_t FieldPower() const{ return fieldProductByIndex.size(); }

    // This returns a string that should be valid Python assuming that the
    // field configuration is given as an array called "fv" and that the
    // Lagrangian parameters are in an array called "lp".
    std::string AsPython() const;

    // This is mainly for debugging.
    std::string AsDebuggingString() const;


  protected:
    double coefficientConstant;
    std::vector< size_t > fieldProductByIndex;
    std::vector< size_t > fieldPowersByIndex;
    std::vector< size_t > parameterIndices;
    double totalCoefficientForFixedScale;

    // This returns doubleToMultiply multiplied by the product of elements of
    // valueVector at the indices given by indexVector.
    double ElementProduct( double doubleToMultiply,
                           std::vector< double > const& valueVector,
                           std::vector< size_t > const& indexVector ) const;
  };






  // This raises the power of the field given by fieldIndex by the number
  // given by powerInt.
  inline void
  ParametersAndFieldsProduct::RaiseFieldPower( size_t const fieldIndex,
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

  // This resets the ParametersAndFieldsProduct to be as if freshly
  // constructed.
  inline void ParametersAndFieldsProduct::ResetValues()
  {
    coefficientConstant = 1.0;
    fieldProductByIndex.clear();
    fieldPowersByIndex.clear();
    parameterIndices.clear();
    totalCoefficientForFixedScale = 1.0;
  }

  // This returns a ParametersAndFieldsProduct that is the partial derivative
  // with respect to the field with index fieldIndex.
  inline ParametersAndFieldsProduct
  ParametersAndFieldsProduct::PartialDerivative(
                                                size_t const fieldIndex ) const
  {
    ParametersAndFieldsProduct returnTerm( *this );
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
    return returnTerm;
  }

  // This returns doubleToMultiply multiplied by the field product using the
  // values in fieldConfiguration.
  inline double
  ParametersAndFieldsProduct::ElementProduct( double doubleToMultiply,
                                      std::vector< double > const& valueVector,
                               std::vector< size_t > const& indexVector ) const
  {
    for( std::vector< size_t >::const_iterator
         elementIndex( indexVector.begin() );
         elementIndex < indexVector.end();
         ++elementIndex )
    {
      doubleToMultiply *= valueVector[ *elementIndex ];
    }
    return doubleToMultiply;
  }

} /* namespace VevaciousPlusPlus */

#endif /* PARAMETERSANDFIELDSPRODUCT_HPP_ */
