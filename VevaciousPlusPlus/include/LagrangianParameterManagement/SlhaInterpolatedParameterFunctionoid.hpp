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

namespace VevaciousPlusPlus
{

  class SlhaInterpolatedParameterFunctionoid
  {
  public:
    SlhaInterpolatedParameterFunctionoid( size_t const indexInValuesVector,
                 LHPC::SLHA::SparseManyIndexedBlock< double > const& slhaBlock,
                                          std::string const& indexString );
    SlhaInterpolatedParameterFunctionoid(
                      SlhaInterpolatedParameterFunctionoid const& copySource );
    virtual ~SlhaInterpolatedParameterFunctionoid();


    size_t IndexInValuesVector() const { return indexInValuesVector; }

    // This returns the value of the functionoid for the given logarithm of the
    // scale.
    double operator()( double const logarithmOfScale ) const
    { return scaleLogarithmPowerCoefficients( logarithmOfScale ); }

    // This re-calculates the coefficients of the polynomial of the logarithm
    // of the scale used in evaluating the functionoid.
    void UpdateForNewSlhaParameters();

    // This is for creating a Python version of the potential.
    virtual std::string PythonParameterEvaluation() const;


  protected:
    LHPC::SLHA::SparseManyIndexedBlock< double > const& slhaBlock;
    std::vector< int > const indexVector;
    SimplePolynomial scaleLogarithmPowerCoefficients;
    size_t indexInValuesVector;
  };





  // This is for creating a Python version of the potential.
  inline std::string
  SlhaInterpolatedParameterFunctionoid::PythonParameterEvaluation() const
  {
    std::vector< double > const&
    scaleCoefficients( scaleLogarithmPowerCoefficients.CoefficientVector() );
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 ) << "( " << scaleCoefficients[ 0 ];
    for( size_t whichPower( 1 );
         whichPower < scaleCoefficients.size();
         ++whichPower )
    {
      stringBuilder << " + ( " << scaleCoefficients[ whichPower ]
      << " ) * lnQ**" << whichPower;
    }
    stringBuilder << " )";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHAINTERPOLATEDPARAMETERFUNCTIONOID_HPP_ */
