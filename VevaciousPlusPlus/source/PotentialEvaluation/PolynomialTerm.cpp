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
      returnStream << " " << *functionoidPointer;
    }
    returnStream
    << " }" << std::endl;
    return std::string( returnStream.str() );
  }

} /* namespace VevaciousPlusPlus */
