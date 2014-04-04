/*
 * PolynomialSum.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  PolynomialSum::PolynomialSum() :
    polynomialTerms()
  {
    // This constructor is just an initialization list.
  }

  PolynomialSum::PolynomialSum( PolynomialSum const& copySource ) :
    polynomialTerms( copySource.polynomialTerms )
  {
    // This constructor is just an initialization list.
  }

  PolynomialSum::~PolynomialSum()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
