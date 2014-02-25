/*
 * PolynomialSum.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALSUM_HPP_
#define POLYNOMIALSUM_HPP_

#include  "../StandardIncludes.hpp"

namespace VevaciousPlusPlus
{
  class PolynomialSum
  {
  public:
    PolynomialSum();
    virtual
    ~PolynomialSum();


    double operator()( std::vector< double > const& fieldConfiguration ) const;


  protected:
  };

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALSUM_HPP_ */
