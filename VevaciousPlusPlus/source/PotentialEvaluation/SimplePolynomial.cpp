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

  SimplePolynomial::SimplePolynomial( size_t const reserveSize,
                                      size_t const leadingPower ) :
    coefficientVector( reserveSize,
                       0.0 ),
    leadingPower( leadingPower )
  {
    // This constructor is just an initialization list.
  }

  SimplePolynomial::SimplePolynomial(
                                std::vector< double > const& coefficientVector,
                                      size_t const leadingPower ) :
    coefficientVector( coefficientVector ),
    leadingPower( leadingPower )
  {
    // This constructor is just an initialization list.
  }

  SimplePolynomial::SimplePolynomial( Eigen::VectorXd const& eigenVector,
                                      size_t const leadingPower,
                                     size_t extraEmptyEntriesAtConstruction ) :
    coefficientVector(),
    leadingPower( leadingPower )
  {
    CopyFromEigen( eigenVector,
                   leadingPower,
                   extraEmptyEntriesAtConstruction );
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
