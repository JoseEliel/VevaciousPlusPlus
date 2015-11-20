/*
 * ParametersAndFieldsProduct.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: bol
 */

#include "BasicFunctions/ParametersAndFieldsProduct.hpp"

namespace VevaciousPlusPlus
{

  ParametersAndFieldsProduct::ParametersAndFieldsProduct() :
    isValid( true ),
    coefficientConstant( 1.0 ),
    fieldProductByIndex(),
    fieldPowersByIndex(),
    parameterIndices(),
    totalCoefficientForFixedScale( 1.0 )
  {
    // This constructor is just an initialization list.
  }

  ParametersAndFieldsProduct::ParametersAndFieldsProduct(
                               ParametersAndFieldsProduct const& copySource ) :
    isValid( copySource.isValid ),
    coefficientConstant( copySource.coefficientConstant ),
    fieldProductByIndex( copySource.fieldProductByIndex ),
    fieldPowersByIndex( copySource.fieldPowersByIndex ),
    parameterIndices( copySource.parameterIndices ),
    totalCoefficientForFixedScale( copySource.totalCoefficientForFixedScale )
  {
    // This constructor is just an initialization list.
  }

  ParametersAndFieldsProduct::~ParametersAndFieldsProduct()
  {
    // This does nothing.
  }


  // This returns a ParametersAndFieldsProduct that is the partial derivative
  // with respect to the field with index fieldIndex.
  ParametersAndFieldsProduct ParametersAndFieldsProduct::PartialDerivative(
                                                size_t const fieldIndex ) const
  {
    if( !(NonZeroDerivative( fieldIndex )) )
    {
      ParametersAndFieldsProduct returnTerm;
      returnTerm.coefficientConstant = 0.0;
      returnTerm.totalCoefficientForFixedScale = 0.0;
      return returnTerm;
    }
    ParametersAndFieldsProduct returnTerm( *this );
    returnTerm.coefficientConstant *= fieldPowersByIndex[ fieldIndex ];
    returnTerm.totalCoefficientForFixedScale
    *= fieldPowersByIndex[ fieldIndex ];
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

  // This returns a string that should be valid Python assuming that the
  // field configuration is given as an array called "fv" and that the
  // Lagrangian parameters are in an array called "lp".
  std::string ParametersAndFieldsProduct::AsPython() const
  {
    if( !isValid )
    {
      return "( 0.0 )";
    }
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 ) << "( " << coefficientConstant;
    for( size_t parameterIndex( 0 );
         parameterIndex < parameterIndices.size();
         ++parameterIndex )
    {
      stringBuilder << " * lp[ " << parameterIndex << " ]";
    }
    for( size_t fieldIndex( 0 );
         fieldIndex < fieldPowersByIndex.size();
         ++fieldIndex )
    {
      if( fieldPowersByIndex[ fieldIndex ] == 1 )
      {
        stringBuilder << " * fv[ " << fieldIndex << " ]";
      }
      if( fieldPowersByIndex[ fieldIndex ] > 1 )
      {
        stringBuilder << " * (fv[ " << fieldIndex << " ])**"
        << fieldPowersByIndex[ fieldIndex ];
      }
    }
    stringBuilder << " )";

    return stringBuilder.str();
  }

  // This is mainly for debugging.
  std::string ParametersAndFieldsProduct::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream
    << "isValid = " << isValid << std::endl
    << "coefficientConstant = " << coefficientConstant << std::endl
    << "fieldProductByIndex = {";
    for( std::vector< size_t >::const_iterator
         fieldIndex( fieldProductByIndex.begin() );
         fieldIndex < fieldProductByIndex.end();
         ++fieldIndex )
    {
      returnStream << " " << *fieldIndex;
    }
    returnStream
    << " }" << std::endl
    << "parameterIndices = {";
    for( std::vector< size_t >::const_iterator
         parameterIndex( parameterIndices.begin() );
         parameterIndex < parameterIndices.end();
         ++parameterIndex )
    {
      returnStream << " " << *parameterIndex;
    }
    returnStream
    << " }" << std::endl
    << "totalCoefficientForFixedScale = " << totalCoefficientForFixedScale
    << std::endl;
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
