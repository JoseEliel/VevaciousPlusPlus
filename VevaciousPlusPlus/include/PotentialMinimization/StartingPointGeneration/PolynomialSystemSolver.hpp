/*
 * PolynomialSystemSolver.hpp
 *
 *  Created on: Nov 20, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALSYSTEMSOLVER_HPP_
#define POLYNOMIALSYSTEMSOLVER_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialSystemSolver
  {
  public:
    typedef std::pair< double, std::vector< size_t > > FactorWithPowers;
    typedef std::vector< FactorWithPowers > PolynomialConstraint;
    typedef std::vector< PolynomialConstraint > ConstraintSystem;

    PolynomialSystemSolver() { ; /* This does nothing. */ }
    virtual ~PolynomialSystemSolver() { ; /* This does nothing. */ }


    // This should solve the system of polynomial constraints and put all
    // the solutions into systemSolutions.
    virtual void operator()( std::vector< PolynomialConstraint > systemToSolve,
             std::vector< std::vector< double > >& systemSolutions ) const = 0;
  };

} /* namespace VevaciousPlusPlus */

#endif /* POLYNOMIALSYSTEMSOLVER_HPP_ */
