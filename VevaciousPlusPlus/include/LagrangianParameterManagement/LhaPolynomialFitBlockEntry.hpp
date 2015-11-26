/*
 * LhaPolynomialFitBlockEntry.hpp
 *
 *  Created on: Oct 28, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LHAPOLYNOMIALFITBLOCKENTRY_HPP_
#define LHAPOLYNOMIALFITBLOCKENTRY_HPP_

#include "CommonIncludes.hpp"
#include "Eigen/Dense"
#include "BasicFunctions/SimplePolynomial.hpp"
#include "LhaInterpolatedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class LhaPolynomialFitBlockEntry : public LhaInterpolatedParameterFunctionoid
  {
  public:
    LhaPolynomialFitBlockEntry( size_t const indexInValuesVector,
                              LHPC::SlhaSimplisticInterpreter const& lhaParser,
                                std::string const& parameterName  );
    LhaPolynomialFitBlockEntry( LhaPolynomialFitBlockEntry const& copySource );
    virtual ~LhaPolynomialFitBlockEntry();


    // This returns the value of the functionoid for the given logarithm of the
    // scale.
    virtual double operator()( double const logarithmOfScale ) const
    { return scaleLogarithmPowerCoefficients( logarithmOfScale ); }

    // This returns the value of the functionoid for the given logarithm of the
    // scale. It ignores the values of the other parameters.
    virtual double operator()( double const logarithmOfScale,
                        std::vector< double > const& interpolatedValues ) const
    { return scaleLogarithmPowerCoefficients( logarithmOfScale ); }

    // This re-calculates the coefficients of the polynomial of the logarithm
    // of the scale used in evaluating the functionoid.
    virtual void UpdateForNewLhaParameters();

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;

    // This is mainly for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    SimplePolynomial scaleLogarithmPowerCoefficients;
  };





  // This is for creating a Python version of the potential.
  inline std::string LhaPolynomialFitBlockEntry::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::vector< double > const&
    scaleCoefficients( scaleLogarithmPowerCoefficients.CoefficientVector() );
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 )
    << PythonIndent( indentationSpaces )
    << "parameterValues[ " << IndexInValuesVector() << " ] = ( "
    << scaleCoefficients[ 0 ];
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

  // This is mainly for debugging.
  inline std::string LhaPolynomialFitBlockEntry::AsDebuggingString() const
  {
    std::stringstream stringBuilder;
    stringBuilder << "IndexInValuesVector() = " << IndexInValuesVector()
    << ", scaleLogarithmPowerCoefficients = "
    << scaleLogarithmPowerCoefficients.AsDebuggingString();
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* LHAPOLYNOMIALFITBLOCKENTRY_HPP_ */
