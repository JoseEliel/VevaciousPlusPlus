/*
 * SlhaInterpolatedParameterFunctionoid.hpp
 *
 *  Created on: Oct 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHAINTERPOLATEDPARAMETERFUNCTIONOID_HPP_
#define SLHAINTERPOLATEDPARAMETERFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "Eigen/Dense"
#include "BasicFunctions/SimplePolynomial.hpp"
#include "SlhaSourcedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaInterpolatedParameterFunctionoid :
                                         public SlhaSourcedParameterFunctionoid
  {
  public:
    SlhaInterpolatedParameterFunctionoid( size_t const indexInValuesVector,
                 LHPC::SLHA::SparseManyIndexedBlock< double > const& slhaBlock,
                                          std::string const& indexString );
    SlhaInterpolatedParameterFunctionoid(
                      SlhaInterpolatedParameterFunctionoid const& copySource );
    virtual ~SlhaInterpolatedParameterFunctionoid();


    // This should return the value of the functionoid for the given logarithm
    // of the scale.
    virtual double operator()( double const logarithmOfScale ) const = 0;

    // This should return the value of the functionoid for the given logarithm
    // of the scale. It should ignore the values of the other parameters.
    virtual double operator()( double const logarithmOfScale,
                   std::vector< double > const& interpolatedValues ) const = 0;

    // This should re-calculate the coefficients of the polynomial of the
    // logarithm of the scale used in evaluating the functionoid.
    virtual void UpdateForNewSlhaParameters() = 0;

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const = 0;


  protected:
    LHPC::SLHA::SparseManyIndexedBlock< double > const& slhaBlock;
    std::vector< int > const indexVector;
  };

} /* namespace VevaciousPlusPlus */

#endif /* SLHAINTERPOLATEDPARAMETERFUNCTIONOID_HPP_ */
