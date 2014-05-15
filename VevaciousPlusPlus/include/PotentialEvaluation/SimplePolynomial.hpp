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
    SimplePolynomial( std::vector< double > const& coefficientVector );
    SimplePolynomial( Eigen::VectorXd const& eigenVector );
    SimplePolynomial( SimplePolynomial const& copySource );
    virtual
    ~SimplePolynomial();


    // This returns the sum of coefficientVector[ p ] * inputValue^p over p.
    double operator()( double const inputValue ) const;

    // This sets coefficientVector to have elements equal to those of
    // eigenVector.
    void CopyFromEigen( Eigen::VectorXd const& eigenVector );

    std::vector< double > const& CoefficientVector() const
    { return coefficientVector; }
    std::vector< double >& CoefficientVector(){ return coefficientVector; }


  protected:
    std::vector< double > coefficientVector;
  };



  // This returns the sum of coefficientVector[ p ] * inputValue^p over p.
  inline double SimplePolynomial::operator()( double const inputValue ) const
  {
    double returnValue( 0.0 );
    for( unsigned int whichPower( 0 );
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
  SimplePolynomial::CopyFromEigen( Eigen::VectorXd const& eigenVector )
  {
    coefficientVector.resize( eigenVector.rows() );
    for( unsigned int whichIndex( 0 );
         whichIndex < coefficientVector.size();
         ++whichIndex )
    {
      coefficientVector[ whichIndex ] = eigenVector( whichIndex );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* SIMPLEPOLYNOMIAL_HPP_ */
