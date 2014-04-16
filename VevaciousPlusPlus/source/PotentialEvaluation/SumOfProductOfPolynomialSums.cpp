/*
 * SumOfProductOfPolynomialSums.cpp
 *
 *  Created on: Apr 16, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  SumOfProductOfPolynomialSums::SumOfProductOfPolynomialSums() :
    sumVector()
  {
    // This constructor is just an initialization list.
  }

  SumOfProductOfPolynomialSums::SumOfProductOfPolynomialSums(
                             SumOfProductOfPolynomialSums const& copySource ) :
    sumVector( copySource.sumVector )
  {
    // This constructor is just an initialization list.
  }

  SumOfProductOfPolynomialSums::~SumOfProductOfPolynomialSums()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
