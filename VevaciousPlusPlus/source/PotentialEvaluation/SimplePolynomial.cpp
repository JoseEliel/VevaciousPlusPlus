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
    coefficientVector(),
    leadingPower( 0 )
  {
    // This constructor is just an initialization list.
  }

  SimplePolynomial::SimplePolynomial( unsigned int const reserveSize,
                                      unsigned int const leadingPower ) :
    coefficientVector( reserveSize,
                       0.0 ),
    leadingPower( leadingPower )
  {
    // This constructor is just an initialization list.
  }

  SimplePolynomial::SimplePolynomial(
                                std::vector< double > const& coefficientVector,
                                      unsigned int const leadingPower ) :
    coefficientVector( coefficientVector ),
    leadingPower( leadingPower )
  {
    // This constructor is just an initialization list.
  }

  SimplePolynomial::SimplePolynomial( Eigen::VectorXd const& eigenVector,
                                      unsigned int const leadingPower ) :
    coefficientVector(),
    leadingPower( leadingPower )
  {
    CopyFromEigen( eigenVector,
                   leadingPower );
  }

  SimplePolynomial::SimplePolynomial( SimplePolynomial const& copySource ) :
    coefficientVector( copySource.coefficientVector ),
    leadingPower( copySource.leadingPower )
  {
    // This constructor is just an initialization list.
  }

  SimplePolynomial::~SimplePolynomial()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
