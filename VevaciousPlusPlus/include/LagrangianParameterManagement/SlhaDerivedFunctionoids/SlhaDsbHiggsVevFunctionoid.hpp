/*
 * SlhaDsbHiggsVevFunctionoid.hpp
 *
 *  Created on: Oct 26, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHADSBHIGGSVEVFUNCTIONOID_HPP_
#define SLHADSBHIGGSVEVFUNCTIONOID_HPP_

#include "../LhaSourcedParameterFunctionoid.hpp"
#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class SlhaDsbHiggsVevFunctionoid : public LhaSourcedParameterFunctionoid
  {
  public:
    SlhaDsbHiggsVevFunctionoid( size_t const indexInValuesVector,
                                size_t const vevIndex,
                                size_t const tanBetaIndex,
                                bool const sinNotCos );
    virtual ~SlhaDsbHiggsVevFunctionoid();


    // This returns the value of the DSB VEV for the neutral component of the
    // functionoid's Higgs doublet: HMIX[3] times either cos(beta) or sin(beta)
    // is returned depending on whether this functionoid is for the doublet
    // which gives mass to the down-type or up-type quarks respectively, with
    // beta calculated from the value of the functionoid evaluating tan(beta).
    virtual double operator()( double const logarithmOfScale,
                       std::vector< double > const& interpolatedValues ) const
    {return ( interpolatedValues[ vevIndex ]
              * (*cosOrSin)( atan( interpolatedValues[ tanBetaIndex ] ) ) ); }

    // This returns the value of the DSB VEV for the neutral component of the
    // functionoid's Higgs doublet: HMIX[3] times either cos(beta) or sin(beta)
    // is returned depending on whether this functionoid is for the doublet
    // which gives mass to the down-type or up-type quarks respectively, with
    // beta calculated from the value of the functionoid evaluating tan(beta).
    double operator()( double const vevEuclideanLength,
                       double const tanBeta ) const
    { return ( vevEuclideanLength * (*cosOrSin)( atan( tanBeta ) ) ); }

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;

    // This is mainly for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    size_t const vevIndex;
    size_t const tanBetaIndex;
    double (*cosOrSin)( double );
    std::string cosOrSinPythonString;
  };





  // This is for creating a Python version of the potential.
  inline std::string SlhaDsbHiggsVevFunctionoid::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 )
    << PythonIndent( indentationSpaces )
    << "parameterValues[ " << IndexInValuesVector()
    << " ] = " << cosOrSinPythonString << "( math.atan( parameterValues[ "
    << tanBetaIndex << " ] ) )";
    return stringBuilder.str();
  }

  // This is mainly for debugging.
  inline std::string
  SlhaDsbHiggsVevFunctionoid::AsDebuggingString() const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << "IndexInValuesVector() = " << IndexInValuesVector() << std::endl;
    stringBuilder << "vevIndex = " << vevIndex << std::endl;
    stringBuilder << "tanBetaIndex = " << tanBetaIndex << std::endl;
    stringBuilder
    << "cosOrSinPythonString = " << cosOrSinPythonString << std::endl;
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHADSBHIGGSVEVFUNCTIONOID_HPP_ */
