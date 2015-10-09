/*
 * PolynomialSum.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BasicFunctions/PolynomialSum.hpp"

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


  // This returns a string that should be valid Python assuming that the
  // field configuration is given as an array called "fv".
  std::string PolynomialSum::AsPython() const
  {
    std::stringstream stringBuilder;
    stringBuilder << "( ";
    if( polynomialTerms.empty() )
    {
      stringBuilder << "0.0";
    }
    else
    {
      for( std::vector< PolynomialTerm >::const_iterator
           polynomialTerm( polynomialTerms.begin() );
           polynomialTerm < polynomialTerms.end();
           ++polynomialTerm )
      {
        if( polynomialTerm != polynomialTerms.begin() )
        {
          stringBuilder << " + ";
        }
        stringBuilder << polynomialTerm->AsPython();
      }
    }
    stringBuilder << " )";

    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
