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
    functionoidProduct()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PolynomialTerm::PolynomialTerm()";
    std::cout << std::endl;/**/
  }

  PolynomialTerm::PolynomialTerm( PolynomialTerm const& copySource ) :
    isValid( copySource.isValid ),
    coefficientConstant( copySource.coefficientConstant ),
    fieldProductByIndex( copySource.fieldProductByIndex ),
    functionoidProduct( copySource.functionoidProduct )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PolynomialTerm::PolynomialTerm( [copy] )";
    std::cout << std::endl;/**/
  }

  PolynomialTerm::~PolynomialTerm()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PolynomialTerm::~PolynomialTerm()";
    std::cout << std::endl;/**/
  }


  // This adds runningParameter to the set of functionoids which multiply
  // coefficientConstant to form the scale-dependent coefficient.
  void
  PolynomialTerm::MultiplyBy( ParameterFunctionoid* const runningParameter,
                              unsigned int const powerInt )
  {
    if( runningParameter == NULL )
    {
      isValid = false;
    }
    else
    {
      for( unsigned int powerCount( 0 );
           powerCount < powerInt;
           ++powerCount )
      {
        functionoidProduct.push_back( runningParameter );
      }
    }
  }
} /* namespace VevaciousPlusPlus */
