/*
 * ProductOfPolynomialSums.cpp
 *
 *  Created on: Apr 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  ProductOfPolynomialSums::ProductOfPolynomialSums() :
    productVector()
  {
    // This constructor is just an initialization list.
  }

  ProductOfPolynomialSums::ProductOfPolynomialSums(
                                  ProductOfPolynomialSums const& copySource ) :
    productVector( copySource.productVector )
  {
    // This constructor is just an initialization list.
  }

  ProductOfPolynomialSums::~ProductOfPolynomialSums()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
