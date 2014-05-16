/*
 * SimplePolynomial.hpp
 *
 *  Created on: May 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SIMPLEPOLYNOMIAL_HPP_
#define SIMPLEPOLYNOMIAL_HPP_

#include "../CommonIncludes.hpp"
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{

  // This is a class to just hold a bunch of coefficients of a polynomial and
  // return the value of the polynomial for a given input.
  class SimplePolynomial
  {
  public:
    SimplePolynomial();
    SimplePolynomial( unsigned int const reserveSize,
                      unsigned int const leadingPower = 0 );
    SimplePolynomial( std::vector< double > const& coefficientVector,
                      unsigned int const leadingPower = 0 );
    SimplePolynomial( Eigen::VectorXd const& eigenVector,
                      unsigned int const leadingPower = 0 );
    SimplePolynomial( SimplePolynomial const& copySource );
    virtual
    ~SimplePolynomial();


    // This returns the sum of coefficientVector[ p ] * inputValue^p over p.
    double operator()( double const inputValue ) const;

    // This sets coefficientVector to have elements equal to those of
    // eigenVector.
    void CopyFromEigen( Eigen::VectorXd const& eigenVector,
                        unsigned int const leadingPower = 0 );

    std::vector< double > const& CoefficientVector() const
    { return coefficientVector; }
    std::vector< double >& CoefficientVector(){ return coefficientVector; }
    unsigned int const& LeadingPower() const{ return leadingPower; }
    unsigned int& LeadingPower(){ return leadingPower; }

  protected:
    std::vector< double > coefficientVector;
    unsigned int leadingPower;
  };



  // This returns the sum of coefficientVector[ p ] * inputValue^p over p.
  inline double SimplePolynomial::operator()( double const inputValue ) const
  {
    double returnValue( 0.0 );
    for( unsigned int whichPower( leadingPower );
         whichPower < coefficientVector.size();
         ++whichPower )
    {
      returnValue += ( coefficientVector[ whichPower ] * pow( inputValue,
                                                              whichPower ) );
    }
    return returnValue;
  }

  // This sets coefficientVector to have elements equal to those of
  // eigenVector.
  inline void
  SimplePolynomial::CopyFromEigen( Eigen::VectorXd const& eigenVector,
                                   unsigned int const leadingPower )
  {
    this->leadingPower = leadingPower;
    coefficientVector.resize( eigenVector.rows() + leadingPower );
    for( unsigned int whichIndex( 0 );
         whichIndex < leadingPower;
         ++whichIndex )
    {
      coefficientVector[ whichIndex ] = 0.0;
    }
    for( unsigned int whichIndex( leadingPower );
         whichIndex < coefficientVector.size();
         ++whichIndex )
    {
      coefficientVector[ whichIndex ] = eigenVector( whichIndex );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* SIMPLEPOLYNOMIAL_HPP_ */
