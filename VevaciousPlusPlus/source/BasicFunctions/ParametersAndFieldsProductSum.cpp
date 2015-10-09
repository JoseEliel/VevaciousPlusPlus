/*
 * ParametersAndFieldsProductSum.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: bol
 */

#include "BasicFunctions/ParametersAndFieldsProductSum.hpp"

namespace VevaciousPlusPlus
{

  ParametersAndFieldsProductSum::ParametersAndFieldsProductSum() :
    parametersAndFieldsProducts()
  {
    // This constructor is just an initialization list.
  }

  ParametersAndFieldsProductSum::ParametersAndFieldsProductSum(
                            ParametersAndFieldsProductSum const& copySource ) :
    parametersAndFieldsProducts( copySource.parametersAndFieldsProducts )
  {
    // This constructor is just an initialization list.
  }

  ParametersAndFieldsProductSum::~ParametersAndFieldsProductSum()
  {
    // This does nothing.
  }


  // This returns a string that should be valid Python assuming that the
  // field configuration is given as an array called "fv" and that the
  // Lagrangian parameters are in an array called "lp".
  std::string ParametersAndFieldsProductSum::AsPython() const
  {
    std::stringstream stringBuilder;
    stringBuilder << "( ";
    if( parametersAndFieldsProducts.empty() )
    {
      stringBuilder << "0.0";
    }
    else
    {
      for( std::vector< PolynomialTerm >::const_iterator
           parametersAndFieldsProduct( parametersAndFieldsProducts.begin() );
           parametersAndFieldsProduct < parametersAndFieldsProducts.end();
           ++parametersAndFieldsProduct )
      {
        if( parametersAndFieldsProduct != parametersAndFieldsProducts.begin() )
        {
          stringBuilder << " + ";
        }
        stringBuilder << parametersAndFieldsProduct->AsPython();
      }
    }
    stringBuilder << " )";

    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
