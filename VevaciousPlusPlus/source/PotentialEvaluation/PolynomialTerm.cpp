/*
 * PolynomialTerm.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  PolynomialTerm::PolynomialTerm( bool const coefficientIsPositive ) :
    isValid( true ),
    coefficientConstant( -1.0 )
  {
    if( coefficientIsPositive )
    {
      coefficientConstant = 1.0;
    }

    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PolynomialTerm::PolynomialTerm( " << coefficientIsPositive << " )";
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

} /* namespace VevaciousPlusPlus */
