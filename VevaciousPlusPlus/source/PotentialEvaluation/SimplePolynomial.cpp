/*
 * SimplePolynomial.cpp
 *
 *  Created on: May 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  SimplePolynomial::SimplePolynomial() :
    coefficientVector()
  {
    // This constructor is just an initialization list.
  }

  SimplePolynomial::SimplePolynomial(
                             std::vector< double > const& coefficientVector ) :
    coefficientVector( coefficientVector )
  {
    // This constructor is just an initialization list.
  }

  SimplePolynomial::SimplePolynomial( Eigen::VectorXd const& eigenVector ) :
    coefficientVector()
  {
    CopyFromEigen( eigenVector );
  }

  SimplePolynomial::SimplePolynomial( SimplePolynomial const& copySource ) :
    coefficientVector( copySource.coefficientVector )
  {
    // This constructor is just an initialization list.
  }

  SimplePolynomial::~SimplePolynomial()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
