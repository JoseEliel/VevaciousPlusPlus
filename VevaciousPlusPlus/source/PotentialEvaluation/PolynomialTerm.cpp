/*
 * PolynomialTerm.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  PolynomialTerm::PolynomialTerm() :
    isValid( true ),
    coefficientConstant( 1.0 ),
    fieldProductByIndex(),
    fieldPowersByIndex(),
    functionoidProduct()
  {
    // This constructor is just an initialization list.
  }

  PolynomialTerm::PolynomialTerm( PolynomialTerm const& copySource ) :
    isValid( copySource.isValid ),
    coefficientConstant( copySource.coefficientConstant ),
    fieldProductByIndex( copySource.fieldProductByIndex ),
    fieldPowersByIndex( copySource.fieldPowersByIndex ),
    functionoidProduct( copySource.functionoidProduct )
  {
    // This constructor is just an initialization list.
  }

  PolynomialTerm::~PolynomialTerm()
  {
    // This does nothing.
  }


  // This prints the polynomial with the current value of coefficientConstant
  // and the powers of the fields with names given by fieldNames to a string
  // and returns it.
  std::string PolynomialTerm::AsStringAtCurrentScale(
                                  std::vector< std::string > const& fieldNames,
                                                   bool const followsOtherTerm,
                                                bool const emptyIfZero ) const
  {
    if( fieldNames.size() < fieldPowersByIndex.size() )
    {
      throw std::out_of_range( "Not enough field names for PolynomialTerm!" );
    }
    std::stringstream returnStream;
    double currentCoefficient( coefficientConstant );
    for( std::vector< ParameterFunctionoid* >::const_iterator
         whichFunctionoid( functionoidProduct.begin() );
         whichFunctionoid < functionoidProduct.end();
         ++whichFunctionoid )
    {
      currentCoefficient *= (*(*whichFunctionoid))();
    }
    if( emptyIfZero
        &&
        ( currentCoefficient == 0.0 ) )
    {
      return std::string( "" );
    }
    if( followsOtherTerm )
    {
      if( currentCoefficient >= 0.0 )
      {
        returnStream << " + ";
      }
      else
      {
        returnStream << " - ";
        currentCoefficient = (-currentCoefficient);
      }
    }
    else if( currentCoefficient < 0.0 )
    {
      returnStream << '-';
      currentCoefficient = (-currentCoefficient);
    }
    returnStream << currentCoefficient;
    for( unsigned int whichField( 0 );
         whichField < fieldPowersByIndex.size();
         ++whichField )
    {
      if( fieldPowersByIndex[ whichField ] > 0 )
      {
        returnStream << " * " << fieldNames[ whichField ] << "^"
        << fieldPowersByIndex[ whichField ];
      }
    }
    return returnStream.str();
  }

  // This is mainly for debugging.
  std::string PolynomialTerm::AsDebuggingString() const
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
      returnStream << " " << *functionoidPointer << " => "
      << (*(*functionoidPointer))();
    }
    returnStream
    << " }" << std::endl;
    return returnStream.str();
  }

  // This returns a string that should be valid Python assuming that the
  // field configuration is given as an array called "fv".
  std::string PolynomialTerm::AsPython() const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 ) << "( " << coefficientConstant;
    for( std::vector< ParameterFunctionoid* >::const_iterator
         runningParameter( functionoidProduct.begin() );
         runningParameter < functionoidProduct.end();
         ++runningParameter )
    {
      stringBuilder << " * " << (*runningParameter)->PythonParameterName();
    }
    for( unsigned int fieldIndex( 0 );
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

} /* namespace VevaciousPlusPlus */
